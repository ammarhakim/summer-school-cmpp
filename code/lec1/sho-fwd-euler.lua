local DataStruct = require "DataStruct"

-- stores particle position and velocity (third component to make
-- plotter happy)
local ptclPosition = DataStruct.DynVector { numComponents = 3 }

omega = 1.0 --frequency
z0 = 0.0 -- initial position
v0 = 1.0 -- initial velocity
tEnd = 2*math.pi -- final time
cflFrac = 0.1 -- time-step factor

-- function to add data to output field
function appendData(t, z, v)
   ptclPosition:appendData(t, {omega^2*z, v, 0.0})
end

-- store initial conditions
appendData(0.0, z0, v0)

-- compute dt
local dt = cflFrac*1.0/omega

-- main loop
local isDone = false
local tCurr = 0.0
while not isDone do
   if (tCurr+dt >= tEnd) then
      isDone = true
      dt = tEnd-tCurr
   end
   -- update solution
   local z1 = z0 + dt*v0
   local v1 = v0 - dt*omega^2*z0

   -- store solution
   appendData(tCurr+dt, z1, v1)

   -- prepare for next time-step
   z0, v0 = z1, v1
   tCurr = tCurr+dt
end

-- write solution
ptclPosition:write("ptclData.bp")
