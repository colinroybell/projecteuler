# Project Euler 566
#
# a,b,c

from math import sqrt

modCache = {}

def lcm(a,b):
    return a * b // gcd(a,b)

def gcd(a,b):
    if a < b:
        a,b = b,a
    r = a % b 
    if (r==0):
        return b
    else:
        return gcd(b,r)    

def xgcd(a, b):
    x0, x1, y0, y1 = 1, 0, 0, 1
    while b != 0:
        q, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return  a, x0, y0        

def oneCase(A,B,C):
    global n,a,b,c,c_int,c_rem,cSquare,w
    print("Running {} {} {}".format(A,B,C))
    n = lcm(A,B)
    cSquare = False
    if (sqrt(C)-int(sqrt(C+0.5))<1e-6):
        cSquare = True
        cSqrt = int(sqrt(C+0.5))
        n = lcm(n,cSqrt)
        c = n//cSqrt
        c_int = c
        c_rem = 0
    else:    
        c = n/sqrt(C)
        c_int = int(c)
        c_rem = c - int(c)

    a = n//A
    b = n//B

    print("{} steps, widths {} {} and {}".format(n,a,b,c))

    w = [n//A, n//B, n/sqrt(C)]

    print(w)

    scheme2()

def to_str(p):
    i,r,orient,s = p
    v = i + c_rem*r
    return "{}+{}e, val={}, or={} s = {}".format(i,r,v,orient,s)

def adv(p):
    i,r,orient,s = p
    v = i + c_rem*r
    if (v>=n):
        v-=n
    if (v<0):
        v+=n    

    rationalCase = False
    if (s<2 or cSquare):
        rationalCase = True

    width = 0  

    flip = False
    if (v > 0 and v < w[s]):
        flip = True
        if (orient == +1):
            width = w[s]-v
        else:
            width = v    
    elif (v == 0 and i == 0 and r ==0 and orient == +1):
        flip = True   
        width = w[s] 
    elif (rationalCase and i == w[s] and r == 0 and orient == -1):
        flip =True
        width = w[s]
    elif (not rationalCase and i == c_int and r == 1 and orient == -1):
        flip = True
        width = w[s]
    else:
        if (orient==+1):
            width = n - v
        else:
            if (v > 0):
                width = v - w[s]     
            else: # special case for 0
                width = n - w[s]       
    if flip:
        i = n - i
        r = - r
        orient = -orient
    else:
        if i == 0:
            i = n
        i -= int(w[s])
        if not rationalCase:
            r -= 1
    return ((i,r,orient,(s+1)%3),width)        

def answer(m):
    it_s = 0
    it_f = []
    for s_cand,f_cand in m.items():
        new_f = []
        s,f = factorise(s_cand,f_cand)
        if it_s == 0:
            it_s, it_f = s,f
        else:    
            g,m1,m2 = xgcd(it_s,s)
            size = it_s * s // (g*g) 
            lcm = it_s * s // g

            s1 = it_s//g
            s2 = s//g

            for f1 in it_f:
                for f2 in f:
                    if (f1 % g == f2 % g):
                        # get into standard form
                        e1 = f1//g
                        e2 = f2//g
                        new_e = ((e1*m2*s2 + e2*m1*s1)%size)*g + (f1 %g)
                        new_f.append(new_e)
            it_s = lcm
            it_f = new_f                
    if (len(it_f)==1):
        return it_s
    else:                  
        it_f_sorted = sorted(it_f)    
        return it_f_sorted[1] # not the initial 0    
            
        
def factoriseCandidate(s,f,n):
    cand = set()
    for i in range(n):
        for j in range(s//n):
            if (i in f != i+j*n in f):
                return None
        if i in f:
            cand.add(i)
    return cand                

def factorise(s,f):
    print("Trying to factorise {} {}".format(s,f))
    k = sorted(f)
    for n in k:
        if (n > 0 and s % n == 0):
            print("Trial {}".format(n))
            candSplit = factoriseCandidate(s,f,n)
            if (candSplit):
                return (n,candSplit)
    return (s,f)            

def scheme2():
    status = [(-1,0)]*n
    modulos = {}
    places = {}
    for start in range(n):
        if status[start] == (-1,0): # not done yet
            candidates = set()
            offsets = set()
            hits = set()
            p_start = (start,0,1,0)
            steps = 0
            p = p_start
            if debug:
              print("\nStarting at {}".format(p))
            candidates.add(0)
            offsets.add(0)
            hits.add(start)
            safeWidth = n
            status[start] = (start,0)
            while (1):
                p,stepWidth = adv(p)
                if (stepWidth < safeWidth):
                    safeWidth = stepWidth

                steps += 1
                if debug:
                  print("{}. Advance to {} with {}".format(steps,to_str(p),safeWidth))
                assert (p not in places)
                places[p] = 1
                  
                if (p == p_start):
#                    print("Done")
                    break
                i,r,orient,s = p
                if (orient == 1):
                    candidates.add(steps)

                if (s== 0 and r == 0 and orient == 1):
                    # Run into another starting position
#                    print("New position at {}".format(steps))
                    assert(status[i]==(-1,0))
                    status[i] = (start,steps)
                    offsets.add(steps)
                    hits.add(i)
#            print("Candidates: {}".format(candidates))
#            print("Offsets: {}".format(offsets))        

#            print("Safe width {}".format(safeWidth))
            # Optimise out the other cases
            if 0:
                for hit in hits:
                    s1, s2 = status[hit]
                    for k in range(1,int(safeWidth+0.999)): # ceiling, but be safe                
                        status[hit+k] = (s1+k,s2)
    #                    print("Optimised out {}".format(hit+k))

            filtered = set()
            for c in candidates:
                blocked = False
                for offset in offsets:
                    if (c+offset)%steps not in candidates:
#                        print("{} blocked by {}".format(c,offset))
                        blocked = True
                        break
                if not blocked:
                    filtered.add(c)        
              
#            print("Filtered: {}".format(filtered))   

        if steps not in modulos:
            modulos[steps] = filtered
        else:
            modulos[steps] &= filtered    
    print(modulos)    
    modString = str(modulos)
    a = answer(modulos)
    print("Answer {}".format(a))
    modCache[modString] = a

def checkSteps(steps):
    for i in range(n):
        p = (i,0,1,0)
        for s in range(steps):
            p,_ = adv(p)
        i,r,orient,s = p
        if (orient < 1):
            print("Fail at {}".format(i))    

debug = 1
oneCase(15,16,17)
checkSteps(2304)
assert(0)
maximum = 40
cases = 0
for A in range(9,maximum-1):
    for B in range(A,maximum):
        for C in range(B,maximum+1):
            oneCase(A,B,C)
            cases+=1

print("{} cases, of which {} are distinct".format(cases,len(modCache)))            