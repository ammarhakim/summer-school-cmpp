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
polyOrder = 1 -- polynomial order (DONT CHANGE THIS: WORKS ONLY FOR P=1)
xvel = 1.0 -- x-direction velocity
yvel = 0.0 -- y-direction velocity
cfl = 0.5/(2+1) -- CFL number
tEnd = 2*math.pi*8
useAntiLimiter = false -- if we should use anti-limiters
rescaleSolution = false -- if we should rescale solution
extraType = "none" -- one of "none", "linear", "exp", "exp0", "patch-fit"
initProfile = "sinWave" -- one of "gaussian", "step", "cylinder", "expTent", "square-hat"

muQuad = 3.0/5.0 -- location of surface quadrature nodes
rMax = 3.0 -- maximum slope/mean-value ratio allowed
cflAL = cfl -- CFL number to use in anti-limiter
singleStepSim = false

----------------------
-- Grids and fields --
----------------------

grid = Grid.RectCart {
   lower = {-math.pi, 0.0},
   upper = {math.pi, 1.0},
   cells = {8, 2},
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

fNodal = DataStruct.Field {
   onGrid = grid,
   numComponents = basis:numBasis(),
   ghost = {1, 1},
}
-- diagnostics
deltaChange = DataStruct.DynVector {  numComponents = 1 }
rescaledCells = DataStruct.DynVector {  numComponents = 1 }
density = DataStruct.DynVector {  numComponents = 1 }
absDist = DataStruct.DynVector {  numComponents = 1 }

--------------
-- Updaters --
--------------

function sinWave(t, xn)
   return math.sin(xn[1])+1
end
function gaussian(t, xn)
   local r2 = (xn[1]-0.5)^2
   return math.exp(-2*xn[1]^2)
end
function cylinder(t, xn)
   local r2 = (xn[1]-0.5)^2 + (xn[2]-0.5)^2
   if r2 < 0.25^2 then
      return 1.0
   end
   return 1.0e-5
end
function step(t, xn)
   local r2 = (xn[1]-0.5)^2
   if r2 < 0.25^2 then
      return 1.0
   end
   return 1.0e-5
end
function squareHat(t, xn)
   local rx2, ry2 = (xn[1]-0.5)^2, (xn[2]-0.5)^2
   if rx2 < 0.25^2 and ry2 < 0.25^2 then
      return 1.0
   end
   return 1.0e-5
end
function expTent(t, xn)
   local r = math.sqrt((xn[1]-0.5)^2 + (xn[2]-0.5)^2)
   return math.exp(-10*r)
end

   
-- select initial profile
initFunc = nil
if initProfile == "gaussian" then
   initFunc = gaussian
elseif initProfile == "cylinder" then
   initFunc = cylinder
elseif initProfile == "step" then
   initFunc = step
elseif initProfile == "expTent" then
   initFunc = expTent
elseif initProfile == "square-hat" then
   initFunc = squareHat
elseif initProfile == "sinWave" then
   initFunc = sinWave
else
   initFunc = nil
end

-- updater to initialize distribution function
initDist = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   evaluate = function (t, xn)
      return initFunc(t, xn)
   end
}

-- volume contribution
function calcVolumeTerm(dx, dy, u, v, f, out)
   out[1] = 0.0 
   out[2] = (3.464101615137754*f[1]*u)/dx 
   out[3] = (3.464101615137754*f[1]*v)/dy 
   out[4] = (3.464101615137754*f[2]*v)/dy+(3.464101615137754*f[3]*u)/dx
end

-- x/y-direction numerical flux
function numericalFluxX(vel, fl, fr, out)
   if vel > 0 then
      local f = fl
      out[1] = 1.224744871391589*f[2]*vel+0.7071067811865475*f[1]*vel 
      out[2] = 1.224744871391589*f[4]*vel+0.7071067811865475*f[3]*vel
   else
      local f = fr
      out[1] = 0.7071067811865475*f[1]*vel-1.224744871391589*f[2]*vel 
      out[2] = 0.7071067811865475*f[3]*vel-1.224744871391589*f[4]*vel
   end
end
function numericalFluxY(vel, fb, ft, out)
   if vel > 0 then
      local f = fb
      out[1] = 1.224744871391589*f[3]*vel+0.7071067811865475*f[1]*vel 
      out[2] = 1.224744871391589*f[4]*vel+0.7071067811865475*f[2]*vel
   else
      local f = ft
      out[1] = 0.7071067811865475*f[1]*vel-1.224744871391589*f[3]*vel
      out[2] = 0.7071067811865475*f[2]*vel-1.224744871391589*f[4]*vel
   end
end

-- evaluate function at specified coordinate
function evalFunc(f, x, y)
   return 1.5*f[4]*x*y+0.8660254037844386*f[3]*y+0.8660254037844386*f[2]*x+0.5*f[1]
end
-- attach floor
function evalFuncFloor(f, x, y)
   return math.max(0.0, evalFunc(f, x, y))
end

-- compute average and slope in X along a specified y constant line
function averageSlopeX(y, f)
   return 0.8660254037844386*f[3]*y+0.5*f[1], 1.5*f[4]*y+0.8660254037844386*f[2]
end

-- compute average and slope in Y along a specified x constant line
function averageSlopeY(x, f)
   return 0.8660254037844386*f[2]*x+0.5*f[1], 1.5*f[4]*x+0.8660254037844386*f[3]
end

-- reconstruct expansion give values and locations
function reconsSurfExpansion(mu1, mu2, v1, v2, out)
   out[1] = -(1.414213562373095*mu2*v1-1.414213562373095*mu1*v2)/(mu1-mu2) 
   out[2] = (1.414213562373095*v1-1.414213562373095*v2)/(1.732050807568877*mu1-1.732050807568877*mu2)
end

function patchFit(f0, f1, x, CFL)
   local r = f1/(f0 + GKYL_EPSILON)
   local val = 0.0
   if x > 0 then
      if r<2.2 then
	 val = math.max(0, math.min(1.0/CFL, math.exp(2*r*x/3)*(1+r*x/3)))
      else
	 val = math.min(1.0/CFL, 6/(3-math.min(2.999, math.abs(r))))
      end
   else
      if r>-2.2 then
	 val = math.max(0, math.min(1.0/CFL, math.exp(2*r*x/3)*(1+r*x/3)))
      else
	 val = math.min(1.0/CFL, 6/(3+math.min(2.999, math.abs(r))))
      end	 
   end
   return val
end

-- anti-limiter function
function limTheta(f0, f1, x, CFL)
   local r = f1/(f0 + GKYL_EPSILON)

   if extraType == "none" then
      val = 1+r*x
   elseif extraType == "linear" then
      val = math.max(0, 1+r*x)
   elseif extraType == "exp" then
      val = math.min(1.0/CFL, math.exp(r*x))
   elseif extraType == "exp0" then
      val = math.max(0, math.min(1.0/CFL, math.exp(2*r*x/3)*(1+r*x/3)))
   elseif extraType == "patch-fit" then
      val = patchFit(f0, f1, x, CFL)
   else
      val = 0.0
   end
   return val
end

-- x/y-direction numerical flux: anti-limiter version
function numericalFluxX_AL(vel, fl, fr, out)
   local mu1, mu2 = -muQuad, muQuad
   
   -- calculate left/right values at quadrature nodes
   local f0, f1 = averageSlopeX(mu1, fl)
   local fq1_L = limTheta(f0, f1, 1, cflAL)*vel*f0

   local f0, f1 = averageSlopeX(mu1, fr)
   local fq1_R = limTheta(f0, f1, -1, cflAL)*vel*f0

   local f0, f1 = averageSlopeX(mu2, fl)
   local fq2_L = limTheta(f0, f1, 1, cflAL)*vel*f0

   local f0, f1 = averageSlopeX(mu2, fr)
   local fq2_R = limTheta(f0, f1, -1, cflAL)*vel*f0

   local fq1, fq2
   if vel > 0 then
      fq1, fq2 = fq1_L, fq2_L
   else
      fq1, fq2 = fq1_R, fq2_R
   end
   reconsSurfExpansion(mu1, mu2, fq1, fq2, out)
end
function numericalFluxY_AL(vel, fb, ft, out) -- SOMETHING BAD IS GOING ON HERE
   local mu1, mu2 = -muQuad, muQuad
   
   -- calculate left/right values at quadrature nodes
   local f0, f1 = averageSlopeY(mu1, fb)
   local fq1_B = limTheta(f0, f1, 1, cflAL)*vel*f0

   local f0, f1 = averageSlopeY(mu1, ft)
   local fq1_T = limTheta(f0, f1, -1, cflAL)*vel*f0

   local f0, f1 = averageSlopeY(mu2, fb)
   local fq2_B = limTheta(f0, f1, 1, cflAL)*vel*f0

   local f0, f1 = averageSlopeY(mu2, ft)
   local fq2_T = limTheta(f0, f1, -1, cflAL)*vel*f0

   local fq1, fq2
   if vel > 0 then
      fq1, fq2 = fq1_B, fq2_B
   else
      fq1, fq2 = fq1_T, fq2_T
   end
   reconsSurfExpansion(mu1, mu2, fq1, fq2, out)
end

-- functions to compute right/left flux contribution
function rightCellSurf(dx, dy, fn, out)
   out[1] = -(1.414213562373095*fn[1])/dx 
   out[2] = (2.449489742783178*fn[1])/dx 
   out[3] = -(1.414213562373095*fn[2])/dx 
   out[4] = (2.449489742783178*fn[2])/dx 
end
function leftCellSurf(dx, dy, fn, out)
   out[1] = (1.414213562373095*fn[1])/dx 
   out[2] = (2.449489742783178*fn[1])/dx 
   out[3] = (1.414213562373095*fn[2])/dx 
   out[4] = (2.449489742783178*fn[2])/dx
end

-- functions to compute top/bottom flux contribution
function topCellSurf(dx, dy, fn, out)
   out[1] = -(1.414213562373095*fn[1])/dy 
   out[2] = -(1.414213562373095*fn[2])/dy 
   out[3] = (2.449489742783178*fn[1])/dy 
   out[4] = (2.449489742783178*fn[2])/dy 
end
function bottomCellSurf(dx, dy, fn, out)
   out[1] = (1.414213562373095*fn[1])/dy 
   out[2] = (1.414213562373095*fn[2])/dy 
   out[3] = (2.449489742783178*fn[1])/dy 
   out[4] = (2.449489742783178*fn[2])/dy 
end

-- vec = vec + fact*rhs
function incrementVec(vec, fact, rhs)
   vec[1] = vec[1] + fact*rhs[1]
   vec[2] = vec[2] + fact*rhs[2]
   vec[3] = vec[3] + fact*rhs[3]
   vec[4] = vec[4] + fact*rhs[4]
end

-- function to take a single forward Euler step
function forwardEuler(dt, fIn, fOut)
   local volTerm = Lin.Vec(4)
   local surfTerm = Lin.Vec(4)
   local fluxCoeff = Lin.Vec(2)
   local dx, dy = grid:dx(1), grid:dx(2)
   local localRange = fIn:localRange()

   local indexer = fIn:genIndexer()
   fOut:copy(fIn)
   -- loop over each cell, accumulating volume contributions
   for idx in localRange:colMajorIter() do
      local fInPtr, fOutPtr = fIn:get(indexer(idx)), fOut:get(indexer(idx))
      calcVolumeTerm(dx, dy, xvel, yvel, fInPtr, volTerm)
      incrementVec(fOutPtr, dt, volTerm)
   end

   local vel = {xvel, yvel} -- velocity vector

   local numericalFluxTbl = {numericalFluxX, numericalFluxY} -- numeric flux function table
   if useAntiLimiter then
      numericalFluxTbl = {numericalFluxX_AL, numericalFluxY_AL} -- numeric flux function table
   end
   local rightCellSurfTbl = {rightCellSurf, topCellSurf} -- functions to compute numerical flux contribution to right cell
   local leftCellSurfTbl = {leftCellSurf, bottomCellSurf} -- functions to compute numerical flux contribution to left cell

   -- accumulate contributions from surface integrals
   for dir = 1, 2 do
      -- lower/upper bounds in direction 'dir': these are edge indices (one more edge than cell)
      local dirLoIdx, dirUpIdx = localRange:lower(dir), localRange:upper(dir)+1
      local perpRange = localRange:shorten(dir) -- range orthogonal to 'dir'

      -- outer loop is over directions orthogonal to 'dir' and inner
      -- loop is over 1D slice in `dir`.
      for idx in perpRange:colMajorIter() do
   	 local idxp, idxm = idx:copy(), idx:copy()

   	 for i = dirLoIdx, dirUpIdx do -- this loop is over edges
   	    idxm[dir], idxp[dir]  = i-1, i -- cell left/right of edge 'i'
   	    local fInL, fInR = fIn:get(indexer(idxm)), fIn:get(indexer(idxp))

   	    numericalFluxTbl[dir](vel[dir], fInL, fInR, fluxCoeff) -- compute numerical flux

   	    local fOutL, fOutR = fOut:get(indexer(idxm)), fOut:get(indexer(idxp))
   	    -- accumulate contribution to right/left cell
   	    rightCellSurfTbl[dir](dx, dy, fluxCoeff, surfTerm)
	    incrementVec(fOutR, -dt, surfTerm)
	    
   	    leftCellSurfTbl[dir](dx, dy, fluxCoeff, surfTerm)
	    incrementVec(fOutL, -dt, surfTerm)
   	 end
      end
   end
end

function applyBc(fIn)
   local localRange = fIn:localRange()
   local indexer = fIn:genIndexer()

   for dir = 1,2 do
      local perpRange = localRange:shorten(dir) -- range orthogonal to 'dir'
      for idx in perpRange:colMajorIter() do
	 local idxGst, idxSkn = idx:copy(), idx:copy()

	 idxSkn[dir] = localRange:lower(dir)
	 idxGst[dir] = localRange:upper(dir)+1
	 local fGst, fSkn = fIn:get(indexer(idxGst)), fIn:get(indexer(idxSkn))
	 fGst[1] = fSkn[1]
	 fGst[2] = fSkn[2]
	 fGst[3] = fSkn[3]
	 fGst[4] = fSkn[4]

	 idxSkn[dir] = localRange:upper(dir)
	 idxGst[dir] = localRange:lower(dir)-1
	 fGst, fSkn = fIn:get(indexer(idxGst)), fIn:get(indexer(idxSkn))
	 fGst[1] = fSkn[1]
	 fGst[2] = fSkn[2]
	 fGst[3] = fSkn[3]
	 fGst[4] = fSkn[4]
      end
   end
end

-- compute values at interior nodes
function calcNodalField(fIn, nodalOut)
   local localRange = fIn:localRange()
   local indexer = fIn:genIndexer()

   local mu1, mu2 = -1/rMax, 1/rMax

   for idx in localRange:colMajorIter() do
      local fInPtr = fIn:get(indexer(idx))
      local out = nodalOut:get(indexer(idx))
      out[1] = evalFunc(fInPtr, mu1, mu1)
      out[2] = evalFunc(fInPtr, mu2, mu1)
      out[3] = evalFunc(fInPtr, mu1, mu2)
      out[4] = evalFunc(fInPtr, mu2, mu2)
   end
end

-- rescale solution to insure positivity at nodes
function rescaleSol(fIn)
   local localRange = fIn:localRange()
   local indexer = fIn:genIndexer()
   local mu1, mu2 = -1/rMax, 1/rMax
   local d2 = 0.0

   local numRescaledCells = 0
   for idx in localRange:colMajorIter() do
      local fInPtr = fIn:get(indexer(idx))
      local f0 = 0.5*fInPtr[1] -- cell average

      local v1, v2, v3, v4 = evalFunc(fInPtr, mu1, mu1),
      evalFunc(fInPtr, mu2, mu1),
      evalFunc(fInPtr, mu1, mu2),
      evalFunc(fInPtr, mu2, mu2)
      local fm = math.min(v1, v2, v3, v4)
      local theta = math.min(1, f0/(f0-fm+GKYL_EPSILON))

      -- compute diagnostics *before* applying change
      local delChange = (1-theta)^2*(fInPtr[2]^2+fInPtr[3]^2+fInPtr[4]^2)
      
      -- modify moments
      fInPtr[2] = theta*fInPtr[2] -- (note no change to cell averages)
      fInPtr[3] = theta*fInPtr[3]
      fInPtr[4] = theta*fInPtr[4]

      if theta < 1 then
	 numRescaledCells = numRescaledCells+1
	 d2 = d2 + delChange
      end
   end
   return numRescaledCells, math.sqrt(d2)
end

-- compute sum abs(distf)
function calcAbsDist(tCurr, fIn)
   local localRange = fIn:localRange()
   local indexer = fIn:genIndexer()

   local fa = 0.0
   for idx in localRange:colMajorIter() do
      local fInPtr = fIn:get(indexer(idx))
      fa = fa+math.abs(0.5*fInPtr[1])
   end
   absDist:appendData(tCurr, {fa})
end

-- compute density
function calcDensity(tCurr, fIn)
   local localRange = fIn:localRange()
   local indexer = fIn:genIndexer()

   local fa = 0.0
   for idx in localRange:colMajorIter() do
      local fInPtr = fIn:get(indexer(idx))
      fa = fa+0.5*fInPtr[1]
   end
   density:appendData(tCurr, {fa})
end

-- apply rescaling limiter
function applyRescaleLimiter(fIn)
   if rescaleSolution then
      return rescaleSol(fIn)
   end
   return 0, 0.0
end

-- take a single time-step with RK3 method
function rk3(tCurr, dt, fIn, fOut)
   local totalRescaledCells, totalDc = 0, 0.0
   -- Stage 1
   forwardEuler(dt, fIn, f1)
   local nrs, dc = applyRescaleLimiter(f1)
   --print(string.format("N = %d; Del = %g", nrs, dc))
   
   totalRescaledCells = totalRescaledCells + nrs
   totalDc = totalDc + dc
   applyBc(f1)

   -- Stage 2
   forwardEuler(dt, f1, fe)
   local nrs, dc = applyRescaleLimiter(fe)
   --print(string.format("N = %d; Del = %g", nrs, dc))
   
   totalRescaledCells = totalRescaledCells + nrs
   totalDc = totalDc + dc   
   f2:combine(3.0/4.0, fIn, 1.0/4.0, fe)
   applyBc(f2)

   -- Stage 3
   forwardEuler(dt, f2, fe)
   local nrs, dc = applyRescaleLimiter(fe)
   --print(string.format("N = %d; Del = %g", nrs, dc))
   
   totalRescaledCells = totalRescaledCells + nrs
   totalDc = totalDc + dc
   fOut:combine(1.0/3.0, fIn, 2.0/3.0, fe)
   applyBc(fOut)

   -- append diagnostic data
   rescaledCells:appendData(tCurr, {totalRescaledCells})
   deltaChange:appendData(tCurr, {totalDc})
   
   return totalRescaledCells, totalDc
end

numApp = 1
-- take a single forward Euler time-step
function rk1(dt, fIn, fOut)
   local nrs  = 0.0
   forwardEuler(dt, fIn, fOut)
   if numApp == 1 then
      local nrs = applyRescaleLimiter(fOut)
      numApp = numApp + 1
   end
   applyBc(fOut)
   return nrs
end

-- run simulation with RK3
function runSimulation(tEnd)
   local vmax = math.abs(xvel)/grid:dx(1) + math.abs(yvel)/grid:dx(2)
   local dt = cfl/vmax
   local tCurr = 0.0
   local step = 1
   local isDone = false
   local totalRescaledCells = 0
   local totalDeltaChange = 0.0   

   local tStart = Time.clock()
   while not isDone do
      if (tCurr+dt >= tEnd) then
	 isDone = true
	 dt = tEnd-tCurr
      end
      print(string.format("Step %d at time %g with dt %g ...", step, tCurr, dt))

      local nrc, dc = rk3(tCurr, dt, f, fNew)
      totalRescaledCells = totalRescaledCells + nrc
      totalDeltaChange = totalDeltaChange + dc

      f:copy(fNew)
       -- compute diagnostics
      calcAbsDist(tCurr+dt, fNew)
      calcDensity(tCurr+dt, fNew)

      step = step+1
      tCurr = tCurr+dt
   end
   tEnd = Time.clock()
   print(string.format("Simulation took %g secs", tEnd-tStart))

   f:write("distf_1.bp", tEnd)
   calcNodalField(f, fNodal)
   fNodal:write("nodal_1.bp", 0.0)
   deltaChange:write("deltaChange_1.bp", 0.0)
   rescaledCells:write("rescaledCells_1.bp", 0.0)
   absDist:write("absDist_1.bp", 0.0)
   density:write("density_1.bp", 0.0)
   
   print("Total number of rescaled cells", totalRescaledCells)
   print("Total delta change", totalDeltaChange)
   print("Averge delta change", totalDeltaChange/step/grid:numCells(1)/grid:numCells(2))
end

-- run simulation with RK1, writing data after each time-step
function singleStep(tEnd)
   local vmax = math.abs(xvel)/grid:dx(1) + math.abs(yvel)/grid:dx(2)
   local dt = cfl/vmax
   local tCurr = 0.0
   local step = 1
   local isDone = false
   local totalRescaledCells = 0

   local tStart = Time.clock()
   while not isDone do
      if (tCurr+dt >= tEnd) then
	 isDone = true
	 dt = tEnd-tCurr
      end
      print(string.format("Step %d at time %g with dt %g ...", step, tCurr, dt))

      totalRescaledCells = totalRescaledCells + rk1(dt, f, fNew)
      f:copy(fNew)

      f:write(string.format("distf_%d.bp", step), tCurr+dt)
      calcNodalField(f, fNodal)
      fNodal:write(string.format("nodal_%d.bp", step), tCurr+dt)

      step = step+1
      tCurr = tCurr+dt
   end
   tEnd = Time.clock()
   print(string.format("Simulation took %g secs", tEnd-tStart))

   print("Total number of rescaled cells", totalRescaledCells)
end

-- initialize simulation
initDist:advance(0.0, {}, {f})
applyBc(f)
rescaleSol(f) -- rescale so solution is positive at control nodes
calcAbsDist(0.0, f)
f:write("distf_0.bp", 0.0)
calcNodalField(f, fNodal)
fNodal:write("nodal_0.bp", 0.0)

-- run simulation
if singleStepSim then
   singleStep(tEnd)
else
   runSimulation(tEnd)
end
