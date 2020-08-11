local DataStruct = require "DataStruct"

-- stores particle position and velocity (third component to make
-- plotter happy)
local ptclPosition = DataStruct.DynVector { numComponents = 3 }
local exactPosition = DataStruct.DynVector { numComponents = 3 }

omega = 1.0 --frequency
z0 = 0.0 -- initial position
v0 = 1.0 -- initial velocity
tEnd = 10*6*math.pi -- final time
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

-- single forward Euler step
function forwardEuler(dt, zold, vold)
   return zold+dt*vold, vold-dt*omega^2*zold
end

-- store initial conditions
appendData(ptclPosition, 0.0, z0, v0)
appendData(exactPosition, 0.0, z0, v0)

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
   
   -- RK stage 1
   zr1, vr1 = forwardEuler(dt, z0, v0)
   -- RK stage 2
   zr2, vr2 = forwardEuler(dt, zr1, vr1)
   zr2 = 3/4*z0 + 1/4*zr2
   vr2 = 3/4*v0 + 1/4*vr2
   -- RK stage 3
   zr3, vr3 = forwardEuler(dt, zr2, vr2)
   z1 = 1/3*z0 + 2/3*zr3
   v1 = 1/3*v0 + 2/3*vr3

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
