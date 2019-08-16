-- Input file to solve transverse magnetic mode Maxwell equations

cfl = 0.75
L = 1.0
-- define some globals constants
lowerx = 0.0
upperx = L
cellsx = 32
kwave = 2
light_speed = 1.0
tperiod = 1/light_speed

-- cell spacing
dx = (upperx-lowerx)/cellsx

grid = Grid.RectCart1D {
   lower = {lowerx},
   upper = {upperx},
   cells = {cellsx},
   periodicDirs = {0},
}

Elc = DataStruct.Field1D {
   -- (electric field is stored on edges)
   onGrid = grid,
   -- [Ex, Ey, Ez]
   numComponents = 3,
   ghost = {1, 1},
}

ElcNew = DataStruct.Field1D {
   onGrid = grid,
   -- [Ex, Ey, Ez]
   numComponents = 3,
   ghost = {1, 1},
}

Mgn = DataStruct.Field1D {
   -- (magnetic field is stored on faces)
   onGrid = grid,
   -- [Bx, By, Bz]
   numComponents = 3,
   ghost = {1, 1},
}

MgnNew = DataStruct.Field1D {
   onGrid = grid,
   -- [Bx, By, Bz]
   numComponents = 3,
   ghost = {1, 1},
}
MgnNew:clear(0.0)

-- create duplicates in case we need to take step again
ElcDup = ElcNew:duplicate()
MgnDup = MgnNew:duplicate()

function initElc(x, y, z)
   local c = light_speed
   local phi = 2*math.pi/L*(kwave*x)
   local knorm = math.sqrt(kwave^2)
   local kxn = kwave/knorm
   local E0 = 1.0/math.sqrt(2.0)
   local Ex = 0.0
   local Ey = E0*math.cos(phi)
   local Ez = 0.0
   local Bx = 0.0
   local By = 0.0
   local Bz = E0/c*math.cos(phi)*kxn
   return Ex, Ey, Ez
end
-- initialize electric field
Elc:set(initElc)
ElcNew:copy(Elc)

function initMgn(x, y, z)
   local c = light_speed
   local phi = 2*math.pi/L*(kwave*x)
   local knorm = math.sqrt(kwave^2)
   local kxn = kwave/knorm
   local E0 = 1.0/math.sqrt(2.0)
   local Ex = 0.0
   local Ey = E0*math.cos(phi)
   local Ez = 0.0
   local Bx = 0.0
   local By = 0.0
   local Bz = E0/c*math.cos(phi)*kxn
   return Bx, By, Bz
end
-- initialize electric field
Mgn:set(initMgn)
MgnNew:copy(Mgn)

-- write initial conditions
Elc:write("electricField_0.h5")
Mgn:write("magneticField_0.h5")

-- updater for electric field
elcUpdate = Updater.EdgeFaceCurl1D {
   onGrid = grid,
   -- ElcNew = Elc + c^2*dt*curl(Mgn)
   alpha = light_speed^2,
   -- speed of propagation
   speed = light_speed,
   -- CFL numbers
   cfl = cfl,
}
-- set input/out arrays (these do not change so set it once)
elcUpdate:setIn( {Elc, Mgn} )
elcUpdate:setOut( {ElcNew} )

-- updater for magnetic field
mgnUpdate = Updater.FaceEdgeCurl1D {
   onGrid = grid,
   -- MgnNew = Mgn - dt*curl(ElcNew)
   alpha = -1.0,
   -- speed of propagation
   speed = light_speed,
   -- CFL numbers
   cfl = cfl,
}
-- set input/out arrays (these do not change so set it once)
mgnUpdate:setIn( {Mgn, ElcNew} )
mgnUpdate:setOut( {MgnNew} )

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)

   local step = 1
   local tCurr = tStart
   local myDt = initDt
   -- main loop
   while true do
      -- copy data in case we need to take step again
      ElcDup:copy(Elc)
      MgnDup:copy(Mgn)

      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
	 myDt = tEnd-tCurr
      end

      Lucee.logInfo (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))

      -- apply BCs on magnetic fields
      Mgn:sync()

      -- update electric field
      elcUpdate:setCurrTime(tCurr)
      elcStatus, elcDtSuggested = elcUpdate:advance(tCurr+myDt)

      -- apply BCs on electric field
      ElcNew:sync()

      -- update magnetic field
      mgnUpdate:setCurrTime(tCurr)
      mgnStatus, mgnDtSuggested = mgnUpdate:advance(tCurr+myDt)

      if ( (elcStatus == false) or (mgnStatus == false) ) then
	 -- time-step too large
	 dtSuggested = math.min(elcDtSuggested, mgnDtSuggested)

	 Lucee.logInfo (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 myDt = dtSuggested
	 -- copy over data from duplicate
	 Elc:copy(ElcDup)
	 Mgn:copy(MgnDup)
      else
	 -- step succeeded, proceed to next step
	 dtSuggested = math.min(elcDtSuggested, mgnDtSuggested)

	 Elc:copy(ElcNew)
	 Mgn:copy(MgnNew)

	 tCurr = tCurr + myDt
	 step = step + 1

	 -- check if done
	 if (tCurr >= tEnd) then
	    break
	 end
      end
   end

   return dtSuggested
end

-- parameters to control time-stepping
tStart = 0.0
tEnd = 8*tperiod

-- compute time-step to move magnetic field to dt/2
dtSuggested = cfl*dx/light_speed
-- advance magnetic field by half time-step
mgnUpdate:setCurrTime(0.0)
mgnStatus, mgnDtSuggested = mgnUpdate:advance(0.5*dtSuggested)
Mgn:copy(MgnNew)

nFrames = 16
tFrame = (tEnd-tStart)/nFrames -- time between frames

tCurr = tStart
for frame = 1, nFrames do
   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   -- advance solution between frames
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   -- write out data
   Elc:write( string.format("electricField_%d.h5", frame), tCurr+tFrame)
   Mgn:write( string.format("magneticField_%d.h5", frame), tCurr+tFrame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end
