import math
#METODO BISECCION

def biseccion(f, a, b, tol=1e-10, maxit=10_000):
   
    fa, fb = f(a), f(b)
    if fa == 0:
        return a, 0, 0
    if fb == 0:
        return b, 0, 0
    if fa*fb > 0:
        raise ValueError("f(a) y f(b) deben tener signos opuestos para bisección.")
    
    it = 0
    while it < maxit:
        m = 0.5*(a+b)
        fm = f(m)
        
        if 0.5*(b-a) < tol:
            return m, it, fm
        
        if fm == 0.0:
            return m, it, fm
        
        if fa*fm < 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
        it += 1
    return 0.5*(a+b), it, f(0.5*(a+b))

#Ejercicio 1

def f1(x):
    return x**3 - 7*x**2 + 14*x - 6

intervalos = {
    "(a) [0, 1]": (0.0, 1.0),
    "(b) [1, 3.2]": (1.0, 3.2),
    "(c) [3.2, 4]": (3.2, 4.0)
}

tol = 1e-10
for nombre, (a,b) in intervalos.items():
    try:
        r, it, fr = biseccion(f1, a, b, tol=tol)
        print(f"{nombre}: raíz ≈ {r:.12f}, f(r) ≈ {fr:.3e}, iteraciones = {it}")
    except ValueError as e:
        print(f"{nombre}: {e}")

#Ejercicio 2


def g2(x):
    return x - 2*math.sin(x)

r2, it2, fr2 = biseccion(g2, 0.0, 2.0, tol=1e-10)
print(f"x ≈ {r2:.12f}, g(x) ≈ {fr2:.3e}, iteraciones = {it2}")

#Ejercicio 3

def h3(x):
    return math.tan(x) - x

r3, it3, fr3 = biseccion(h3, 4.4, 4.6, tol=1e-10)
print(f"x ≈ {r3:.12f}, h(x) ≈ {fr3:.3e}, iteraciones = {it3}")

#Ejercicio 4

print(f"f(-2) = {q4(-2)}")
print(f"f(0)  = {q4(0)}")
prod = q4(-2)*q4(0)
print("Producto f(a)*f(b) =", prod, "->", "válido" if prod<0 else "inválido")


#Ejercicio 5

def f5(x):
    return (x+3)*(x+1)*x*(x-1)*(x-3)

intervalos5 = {
    "(a) [-1.5, 2.5]": (-1.5, 2.5),
    "(b) [-0.5, 2.4]": (-0.5, 2.4),
    "(c) [-0.5, 3]":   (-0.5, 3.0),
    "(d) [-3, -0.5]":  (-3.0, -0.5),
}
zeros = [-3, -1, 0, 1, 3]

def mas_cercano(x, xs):
    return min(xs, key=lambda z: abs(x - z))

for nombre, (a,b) in intervalos5.items():
    try:
        r, it, fr = biseccion(f5, a, b, tol=1e-10)
        z = mas_cercano(r, zeros)
        print(f"{nombre}: converge a ≈ {r:.12f} (cero más cercano: {z}), f(r) ≈ {fr:.1e}")
    except ValueError as e:
        print(f"{nombre}: {e}")

#Ejercicio Aplicado 1


def area_segmento(h, r):
    if h <= 0: 
        return 0.0
    if h >= r:
        return 0.5*math.pi*r*r  
    return r*r*math.acos((r-h)/r) - (r-h)*math.sqrt(2*r*h - h*h)

def profundidad_por_biseccion(L, r, V, tol=1e-2):
    def F(h):
        return L*area_segmento(h, r) - V
    return biseccion(F, 0.0, r, tol=tol)


L = 10.0  
r = 1.0   
V = 12.4  

h, it_h, Fh = profundidad_por_biseccion(L, r, V, tol=1e-2)
print(f"h ≈ {h:.3f} cm, residuo ≈ {Fh:.3e}, iteraciones = {it_h}")

#Ejercicio Aplicado 2

def s_t(t, s0, m, k, g=9.81):
    return s0 - (m*g/k)*t + (m*m*g)/(k*k) * (1 - math.exp(-(k/m)*t))

def tiempo_impacto(s0, m, k, g=9.81, tol=1e-2):
   
    T = 1.0
    while s_t(T, s0, m, k, g) > 0 and T < 1e6:
        T *= 2.0
    if T >= 1e6:
        raise RuntimeError("No se encontró un intervalo adecuado.")
    f = lambda t: s_t(t, s0, m, k, g)
    r, it, fr = biseccion(f, 0.0, T, tol=tol)
    return r, it, fr, T

s0 = 300.0  
m  = 0.25   
k  = 0.1    

t_hit, it_hit, s_res, T_scan = tiempo_impacto(s0, m, k, tol=1e-2)
print(f"t_impacto ≈ {t_hit:.3f} s, s(t) ≈ {s_res:.3e}, iteraciones = {it_hit}, T_explorado={T_scan}")


#Ejercicio TEORICO 1

a, b = 1.0, 2.0
tol = 1e-10
cota = math.ceil(math.log2((b-a)/tol))
print(f"Cota mínima de iteraciones: n ≥ {cota}")

def f_et1(x):
    return x**3 - x - 1

r_et1, it_et1, fr_et1 = biseccion(f_et1, a, b, tol=tol)
print(f"Raíz ≈ {r_et1:.12f}, f(r) ≈ {fr_et1:.3e}, iteraciones reales = {it_et1}")