import Jlkep as jk

r = [-400, 6599.206,300]
v = [15.045803, 1, 5]

rf, vf = jk.propagate_lagrangian(r,v,100,jk.MU_EARTH)
print(rf)
