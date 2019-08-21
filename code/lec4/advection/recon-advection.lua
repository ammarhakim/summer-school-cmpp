-- Gkyl ------------------------------------------------------------------------
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Grid = require "Grid"
local DataStruct = require "DataStruct"
local Time = require "Lib.Time"
local Lin = require "Lib.Linalg"
local Basis = require "Basis"
local Updater = require "Updater"

-- Simulation parameters
polyOrder = 2 -- polynomial order
cfl = 0.5/(2*polyOrder+1) -- CFL number
nCell = 8
tEnd = 2*math.pi*100

----------------------
-- Grids and fields --
----------------------

local grid = Grid.RectCart {
   lower = {-math.pi},
   upper = {math.pi},
   cells = {nCell},
   periodicDirs = {1},
}
-- basis functions
local basis = Basis.CartModalSerendipity { ndim = grid:ndim(), polyOrder = polyOrder }

-- fields
local f = DataStruct.Field {
   onGrid = grid,
   numComponents = basis:numBasis(),
   ghost = {1, 1},
}
local fe = DataStruct.Field {
   onGrid = grid,
   numComponents = basis:numBasis(),
   ghost = {1, 1},
}
local f1 = DataStruct.Field {
   onGrid = grid,
   numComponents = basis:numBasis(),
   ghost = {1, 1},
}
local f2 = DataStruct.Field {
   onGrid = grid,
   numComponents = basis:numBasis(),
   ghost = {1, 1},
}
local fNew = DataStruct.Field {
   onGrid = grid,
   numComponents = basis:numBasis(),
   ghost = {1, 1},
}

-- to store integrate density
local density = DataStruct.DynVector { numComponents = 1 }
local fSquare = DataStruct.DynVector { numComponents = 1 }

--------------
-- Updaters --
--------------

local function applyBc(fld)
   fld:sync()
end

-- projection to apply ICs
local initDist = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   evaluate = function (t, xn)
      return math.sin(xn[1])
      --return math.exp(-2*xn[1]^2)
   end
}
initDist:advance(0.0, {}, {f})
applyBc(f)
f:write("f_0.bp", 0.0)

-- integrated density
local densityCalc = Updater.CartFieldIntegratedQuantCalc {
   onGrid = grid,
   basis = basis,
   numComponents = 1,
   quantity = "V"
}
local fSquareCalc = Updater.CartFieldIntegratedQuantCalc {
   onGrid = grid,
   basis = basis,
   numComponents = 1,
   quantity = "V2"
}

local function calcDiag(tm, fld)
   densityCalc:advance(tm, { fld }, { density })
   fSquareCalc:advance(tm, { fld }, { fSquare })
end
calcDiag(0.0, f)
density:write("density_0.bp", 0.0)
fSquare:write("fSquare_0.bp", 0.0)

-- function encoding stencil
local stencilFunc = {}
stencilFunc[0] = function (dt, dx, fOut, fL, f1, fR)
   fOut[1] = f1[1]-((fR[1]-fL[1])*dt)/(2*dx) 
end

stencilFunc[1] = function (dt, dx, fOut, fL, f1, fR)
   fOut[1] = f1[1]-dt*((-(0.5773502691896258*fR[2])/dx)-(0.5773502691896258*fL[2])/dx+(1.154700538379252*f1[2])/dx+(0.5*fR[1])/dx-(0.5*fL[1])/dx) 
   fOut[2] = f1[2]-dt*((-(1.0*fR[2])/dx)+fL[2]/dx+(0.8660254037844386*fR[1])/dx+(0.8660254037844386*fL[1])/dx-(1.732050807568877*f1[1])/dx)
end

stencilFunc[2] = function (dt, dx, fOut, fL, f1, fR)
   fOut[1] = f1[1]-dt*((0.489139870078079*fR[3])/dx-(0.489139870078079*fL[3])/dx-(0.7036456405748563*fR[2])/dx-(0.7036456405748563*fL[2])/dx+(1.407291281149713*f1[2])/dx+(0.5*fR[1])/dx-(0.5*fL[1])/dx) 
   fOut[2] = f1[2]-dt*((0.8472151069828725*fR[3])/dx+(0.8472151069828725*fL[3])/dx+(1.694430213965745*f1[3])/dx-(1.21875*fR[2])/dx+(1.21875*fL[2])/dx+(0.8660254037844386*fR[1])/dx+(0.8660254037844386*fL[1])/dx-(1.732050807568877*f1[1])/dx) 
   fOut[3] = f1[3]-dt*((1.09375*fR[3])/dx-(1.09375*fL[3])/dx-(1.573399484396763*fR[2])/dx-(1.573399484396763*fL[2])/dx-(4.599167723621307*f1[2])/dx+(1.118033988749895*fR[1])/dx-(1.118033988749895*fL[1])/dx)
end

stencilFunc[3] = function (dt, dx, fOut, fL, f1, fR)
   fOut[1] = f1[1]-dt*((-(0.3779644730092272*fR[4])/dx)-(0.3779644730092272*fL[4])/dx+(0.7559289460184544*f1[4])/dx+(0.6987712429686844*fR[3])/dx-(0.6987712429686844*fL[3])/dx-(0.7577722283113838*fR[2])/dx-(0.7577722283113838*fL[2])/dx+(1.515544456622768*f1[2])/dx+(0.5*fR[1])/dx-(0.5*fL[1])/dx) 
   fOut[2] = f1[2]-dt*((-(0.6546536707079771*fR[4])/dx)+(0.6546536707079771*fL[4])/dx+(1.210307295689818*fR[3])/dx+(1.210307295689818*fL[3])/dx+(2.420614591379636*f1[3])/dx-(1.3125*fR[2])/dx+(1.3125*fL[2])/dx+(0.8660254037844386*fR[1])/dx+(0.8660254037844386*fL[1])/dx-(1.732050807568877*f1[1])/dx) 
   fOut[3] = f1[3]-dt*((-(0.8451542547285166*fR[4])/dx)-(0.8451542547285166*fL[4])/dx+(1.690308509457033*f1[4])/dx+(1.5625*fR[3])/dx-(1.5625*fL[3])/dx-(1.694430213965745*fR[2])/dx-(1.694430213965745*fL[2])/dx-(4.357106264483343*f1[2])/dx+(1.118033988749895*fR[1])/dx-(1.118033988749895*fL[1])/dx) 
   fOut[4] = f1[4]-dt*((-(1.0*fR[4])/dx)+fL[4]/dx+(1.84877493221863*fR[3])/dx+(1.84877493221863*fL[3])/dx-(8.134609701761972*f1[3])/dx-(2.00487686654318*fR[2])/dx+(2.00487686654318*fL[2])/dx+(1.322875655532295*fR[1])/dx+(1.322875655532295*fL[1])/dx-(2.645751311064591*f1[1])/dx)
end

stencilFunc[4] = function (dt, dx, fOut, fL, f1, fR)
   fOut[1] = f1[1]-dt*((0.279296875*fR[5])/dx-(0.279296875*fL[5])/dx-(0.617883327946725*fR[4])/dx-(0.617883327946725*fL[4])/dx+(1.23576665589345*f1[4])/dx+(0.8123215699510955*fR[3])/dx-(0.8123215699510955*fL[3])/dx-(0.7870907966686697*fR[2])/dx-(0.7870907966686697*fL[2])/dx+(1.574181593337339*f1[2])/dx+(0.5*fR[1])/dx-(0.5*fL[1])/dx) 
   fOut[2] = f1[2]-dt*((0.4837563778952138*fR[5])/dx+(0.4837563778952138*fL[5])/dx+(0.9675127557904275*f1[5])/dx-(1.07020531715347*fR[4])/dx+(1.07020531715347*fL[4])/dx+(1.406982231239413*fR[3])/dx+(1.406982231239413*fL[3])/dx+(2.813964462478826*f1[3])/dx-(1.36328125*fR[2])/dx+(1.36328125*fL[2])/dx+(0.8660254037844386*fR[1])/dx+(0.8660254037844386*fL[1])/dx-(1.732050807568877*f1[1])/dx) 
   fOut[3] = f1[3]-dt*((0.6245267984032616*fR[5])/dx-(0.6245267984032616*fL[5])/dx-(1.381629123452673*fR[4])/dx-(1.381629123452673*fL[4])/dx+(2.763258246905345*f1[4])/dx+(1.81640625*fR[3])/dx-(1.81640625*fL[3])/dx-(1.75998852581561*fR[2])/dx-(1.75998852581561*fL[2])/dx-(4.225989640783614*f1[2])/dx+(1.118033988749895*fR[1])/dx-(1.118033988749895*fL[1])/dx) 
   fOut[4] = f1[4]-dt*((0.7389500732074931*fR[5])/dx+(0.7389500732074931*fL[5])/dx+(1.477900146414986*f1[5])/dx-(1.634765625*fR[4])/dx+(1.634765625*fL[4])/dx+(2.149200858704158*fR[3])/dx+(2.149200858704158*fL[3])/dx-(7.533757848790918*f1[3])/dx-(2.082446507213006*fR[2])/dx+(2.082446507213006*fL[2])/dx+(1.322875655532295*fR[1])/dx+(1.322875655532295*fL[1])/dx-(2.645751311064591*f1[1])/dx) 
   fOut[5] = f1[5]-dt*((0.837890625*fR[5])/dx-(0.837890625*fL[5])/dx-(1.853649983840175*fR[4])/dx-(1.853649983840175*fL[4])/dx-(12.16720789870719*f1[4])/dx+(2.436964709853287*fR[3])/dx-(2.436964709853287*fL[3])/dx-(2.361272390006008*fR[2])/dx-(2.361272390006008*fL[2])/dx-(5.669760065401246*f1[2])/dx+(1.5*fR[1])/dx-(1.5*fL[1])/dx)
end

-- take a single forward Euler step
local function forwardEuler(dt, fIn, fOut)
   local dx = grid:dx(1)

   local localRange = fIn:localRange()
   local indexer = fIn:indexer()
   
   -- loop over each cell, accumulating contributions
   for i = localRange:lower(1), localRange:upper(1) do
      local fInPtr1 = fIn:get(indexer(i))
      local fInPtrL, fInPtrR = fIn:get(indexer(i-1)), fIn:get(indexer(i+1))
      local fOutPtr = fOut:get(indexer(i))
      -- just call stencil function to update solution
      stencilFunc[polyOrder](dt, dx, fOutPtr, fInPtrL, fInPtr1, fInPtrR)
   end
end

-- take a single time-step with RK3 method
local function rk3(tCurr, dt, fIn, fOut)
   -- Stage 1
   forwardEuler(dt, fIn, f1)
   applyBc(f1)

   -- Stage 2
   forwardEuler(dt, f1, fe)
   f2:combine(3.0/4.0, fIn, 1.0/4.0, fe)
   applyBc(f2)

   -- Stage 3
   forwardEuler(dt, f2, fe)
   fOut:combine(1.0/3.0, fIn, 2.0/3.0, fe)
   applyBc(fOut)
end

-- run simulation with RK3
local function runSimulation(tEnd)
   local tCurr = 0.0
   local step = 1
   local dt = cfl*grid:dx(1)
   
   while not isDone do
      if (tCurr+dt >= tEnd) then
	 isDone = true
	 dt = tEnd-tCurr
      end
      print(string.format("Step %d at time %g with dt %g ...", step, tCurr, dt))
      rk3(tCurr, dt, f, fNew)
      f:copy(fNew)
      calcDiag(tCurr+dt, f) -- compute diagnostics
      step = step+1
      tCurr = tCurr+dt
   end
   f:write("f_1.bp", tEnd)
   density:write("density_1.bp", tEnd)
   fSquare:write("fSquare_1.bp", tEnd)
end

runSimulation(tEnd)
