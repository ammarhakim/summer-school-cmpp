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
nCell = 32
tEnd = 2*math.pi

----------------------
-- Grids and fields --
----------------------

grid = Grid.RectCart {
   lower = {-math.pi},
   upper = {math.pi},
   cells = {nCell},
   periodicDirs = {1},
}
-- basis functions
basis = Basis.CartModalSerendipity { ndim = grid:ndim(), polyOrder = polyOrder }

-- fields
f = DataStruct.Field {
   onGrid = grid,
   numComponents = basis:numBasis(),
   ghost = {1, 1},
}
fe = DataStruct.Field {
   onGrid = grid,
   numComponents = basis:numBasis(),
   ghost = {1, 1},
}
f1 = DataStruct.Field {
   onGrid = grid,
   numComponents = basis:numBasis(),
   ghost = {1, 1},
}
f2 = DataStruct.Field {
   onGrid = grid,
   numComponents = basis:numBasis(),
   ghost = {1, 1},
}
fNew = DataStruct.Field {
   onGrid = grid,
   numComponents = basis:numBasis(),
   ghost = {1, 1},
}

-- to store integrate density
density = DataStruct.DynVector { numComponents = 1 }
fSquare = DataStruct.DynVector { numComponents = 1 }

--------------
-- Updaters --
--------------

function applyBc(fld)
   fld:sync()
end

-- projection to apply ICs
initDist = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   evaluate = function (t, xn)
      if math.abs(xn[1]) < math.pi/4 then
	 return 1.0
      end
      return 0.5
   end
}
initDist:advance(0.0, {}, {f})
f:write("f_0.bp", 0.0)
f:sync()

-- integrated density
densityCalc = Updater.CartFieldIntegratedQuantCalc {
   onGrid = grid,
   basis = basis,
   numComponents = 1,
   quantity = "V"
}
fSquareCalc = Updater.CartFieldIntegratedQuantCalc {
   onGrid = grid,
   basis = basis,
   numComponents = 1,
   quantity = "V2"
}

function calcDiag(tm, fld)
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

-- take a single forward Euler step
function forwardEuler(dt, fIn, fOut)
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
function rk3(tCurr, dt, fIn, fOut)
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
function runSimulation(tEnd)
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
