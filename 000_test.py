import python_wrapper as pw


sigmax = 0.2;
sigmay = 0.25;
x=1.
y=20


part = type('', (), {})() #create empty object
part.x = x
part.y = y
part.px = 0
part.py = 0

pw.weak_strong_single_particle(part, sigmax, sigmay)

print "Ex=%.2e Ey=%.2e\n"%(part.px, part.py)
