local DataStruct = require "DataStruct"

-- stores particle position and velocity (third component to make
-- plotter happy)
local ptclPosition = DataStruct.DynVector { numComponents = 3 }
local exactPosition = DataStruct.DynVector { numComponents = 3 }

omega = 1.0 --frequency
z0 = 0.0 -- initial position
v0 = 1.0 -- initial velocity
tEnd = 10*2*math.pi -- final time
cflFrac = 0.25 -- time-step factor

-- function to add data to output field
function appendData(fld, t, z, v)
   fld:appendData(t, {omega^2*z, v, 0.0})
end

local a, b = z0, v0/omega
-- function to compute exact solution
function exactSolution(t)
   local z = a*math.cos(omega*t) + b*math.sin(omega*t)
   local v = -omega*a*math.sin(omega*t) + omega*b*math.cos(omega*t)
   return z, v
end

-- store initial conditions
appendData(ptclPosition, 0.0, z0, v0)
appendData(exactPosition, 0.0, z0, v0)

-- compute dt
dt = cflFrac*1.0/omega
print(string.format("Using dt*omega %g ...", dt*omega))

-- main loop
isDone = false
tCurr = 0.0
while not isDone do
   if (tCurr+dt >= tEnd) then
      isDone = true
      dt = tEnd-tCurr
   end
   
   -- update solution
   local jfact = 1/(1+omega^2*dt^2/4)   
   local z1 = jfact*((1-omega^2*dt^2/4)*z0 + dt*v0)
   local v1 = jfact*(-omega^2*dt*z0 + (1-omega^2*dt^2/4)*v0)

   -- store solution ..
   appendData(ptclPosition, tCurr+dt, z1, v1)
   -- .. and also exact solution
   local zEx, vEx = exactSolution(tCurr+dt)
   appendData(exactPosition, tCurr+dt, zEx, vEx)

   -- prepare for next time-step
   z0, v0 = z1, v1
   tCurr = tCurr+dt
end

-- write solution
ptclPosition:write("ptclData.bp")
exactPosition:write("exactData.bp")
