import Jlkep as jk

r = [0, 6599.206,0]
v = [9.045803, 0, 0]

rf, vf = jk.propagate_lagrangian(r,v,100,jk.MU_EARTH)
print(rf)
