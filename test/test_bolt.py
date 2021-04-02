from libdenavit.connections import Bolt

b = Bolt('3/4','GroupA-N',hole_type='OVS',surface_type='ClassA')

print(b.d)
print(b.dh)
print(b.Ab)
print(b.rn_bolt_shear(2))
print(b.rn_slip(2))

