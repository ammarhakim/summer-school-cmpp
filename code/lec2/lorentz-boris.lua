local DataStruct = require "DataStruct"
local math = require("math")
local ffi = require("ffi")

charge = 1.0 -- particle charge
mass = 1.0 -- particle mass
Bmax = 1.0 -- estimate for maximum B field
x0, y0, z0 = 1.0, 0.0, 0.0 -- initial position
vx0, vy0, vz0 = 0.0, 0.0, 0.0 -- initial velocity
tEnd = 2*math.pi -- end-time for simulation
cflFrac = 0.1 -- time-step fraction of CFL

-- function to compute electromagentic field
-- returns Ex, Ey, Ez, Bx, By, Bz
function emField(t, x, y, z)
   return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
end

-- Simulation code: no need to modify

ffi.cdef [[
typedef struct { double x[4], v[4]; } particle_t;
]]

function borisPush(t, dt, q, m, ptcl)
   local qmdt = 0.5*q/m*dt
   local x, y, z = ptcl.x[0], ptcl.x[1], ptcl.x[2]
   local E_0, E_1, E_2, B_0, B_1, B_2 = emField(t, x, y, z)

   -- compute t and s vectors
   local t_0 = qmdt*B_0
   local t_1 = qmdt*B_1
   local t_2 = qmdt*B_2

   local tNorm = t_0*t_0 + t_1*t_1 + t_2*t_2
   local s_0 = 2*t_0/(1+tNorm)
   local s_1 = 2*t_1/(1+tNorm)
   local s_2 = 2*t_2/(1+tNorm)

   -- half-step electric field update
   local vm_0 = ptcl.v[0] + qmdt*E_0
   local vm_1 = ptcl.v[1] + qmdt*E_1
   local vm_2 = ptcl.v[2] + qmdt*E_2

   -- rotation around magnetic field
   -- (first compute cp = vm X t)
   local cp_0 =  vm_1*t_2 - vm_2*t_1;
   local cp_1 =  vm_2*t_0 - vm_0*t_2;
   local cp_2 =  vm_0*t_1 - vm_1*t_0;

   local vprime_0 = vm_0 + cp_0
   local vprime_1 = vm_1 + cp_1
   local vprime_2 = vm_2 + cp_2
      
   -- (first compute cp = vprime X s)
   cp_0 =  vprime_1*s_2 - vprime_2*s_1;
   cp_1 =  vprime_2*s_0 - vprime_0*s_2;
   cp_2 =  vprime_0*s_1 - vprime_1*s_0;

   local vp_0 = vm_0 + cp_0
   local vp_1 = vm_1 + cp_1
   local vp_2 = vm_2 + cp_2

   -- half-step electric field update: this gives final particle velocity
   ptcl.v[0] = vp_0 + qmdt*E_0
   ptcl.v[1] = vp_1 + qmdt*E_1
   ptcl.v[2] = vp_2 + qmdt*E_2

   -- update particle position
   ptcl.x[0] = ptcl.x[0] + dt*ptcl.v[0]
   ptcl.x[1] = ptcl.x[1] + dt*ptcl.v[1]
   ptcl.x[2] = ptcl.x[2] + dt*ptcl.v[2]
end

-- function to add data to output field
function appendData(fld, t, z1, z2, z3)
   fld:appendData(t, {z1, z2, z3})
end

-- stores particle position and velocity
local ptclPosition = DataStruct.DynVector { numComponents = 3 }
local ptclVelocity = DataStruct.DynVector { numComponents = 3 }

-- compute dt
omega0 = math.abs(charge)*Bmax/mass
dt = cflFrac*0.5/omega0
-- particle structure
ptcl = ffi.new("particle_t")

ptcl.x[0], ptcl.x[1], ptcl.x[2] = x0, y0, z0
ptcl.v[0], ptcl.v[1], ptcl.v[2] = vx0, vy0, vz0

appendData(ptclPosition, 0.0, ptcl.x[0], ptcl.x[1], ptcl.x[2])
appendData(ptclVelocity, 0.0, ptcl.v[0], ptcl.v[1], ptcl.v[2])

-- main loop
isDone = false
tCurr = 0.0
while not isDone do
   if (tCurr+dt >= tEnd) then
      isDone = true
      dt = tEnd-tCurr
   end
   
   -- update solution
   borisPush(tCurr, dt, charge, mass, ptcl)

   -- store solution
   appendData(ptclPosition, tCurr+dt, ptcl.x[0], ptcl.x[1], ptcl.x[2])
   appendData(ptclVelocity, tCurr+dt, ptcl.v[0], ptcl.v[1], ptcl.v[2])

   tCurr = tCurr+dt
end

ptclPosition:write("ptclPosition.bp")
ptclVelocity:write("ptclVelocity.bp")
