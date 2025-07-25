import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm

#Constante fisica
epsilon0 = 8.854e-12  #Permitividad del vacio





# -----------------------
#FUNCIONES DE CONSTANTE USO
# -----------------------
#Funcion para generar cargas ficticias y puntos de contorno de un conductor circular
def generar_cargas_y_control(n_cargas, radioConductor, centro, radiocargasficticias):
    xc, yc = centro
    ang = np.linspace(0, 2 * np.pi, n_cargas, endpoint=False)
    
    #Puntos de carga dentro del conductor
    r_carga = radiocargasficticias * radioConductor
    cargas = np.stack([xc + r_carga * np.cos(ang), yc + r_carga * np.sin(ang)], axis=1)
    
    #Puntos de contorno sobre el borde del conductor
    r_contorno = radioConductor
    control = np.stack([xc + r_contorno * np.cos(ang), yc + r_contorno * np.sin(ang)], axis=1)
    
    return cargas, control

#Valores del "coeficiente del potencial": P[i,j] = 1 / |ri - rj|
def construir_matriz_P(todas_cargas, todos_control):
    N = len(todos_control)
    P = np.zeros((N, N))
    for i in range(N):
        ri = todos_control[i]
        for j in range(N):
            rj = todas_cargas[j]
            dist = np.linalg.norm(ri - rj)
            P[i, j] = 1 / (4 * np.pi * epsilon0 * dist)
    return P

#Calcular el potencial escalar electrico en una malla de puntos
def calcular_potencial_malla(todas_cargas, Q, X, Y):
    V = np.zeros_like(X)
    for q, rj in zip(Q, todas_cargas):
        R = np.sqrt((X - rj[0])**2 + (Y - rj[1])**2)
        V += q / (4 * np.pi * epsilon0 * R)
    return V
# -----------------------
#FUNCIONES DE CONSTANTE USO
# -----------------------





# -----------------------
#CONFIGURACION PRINCIPAL
# -----------------------
#Valores fijos
cargas_ficticias_por_conductor = 22
radio_cargas_ficticias = 0.7

#Valores variables
cantidad_de_conductores = 1
voltaje_fase = 500000 #Voltios
radio_conductores = 0.014 #Metros

#Posiciones (se debe tener la misma cantidad de posiciones que de la variable cantidad_de_conductores)
centros_conductores = [(0, 0)]
potenciales = [voltaje_fase] * cantidad_de_conductores  #Potencial escalar electrico en cada conductor
# -----------------------
#CONFIGURACION PRINCIPAL
# -----------------------





# -----------------------
#METODO DE SIMULACION DE CARGAS
# -----------------------
#Generar cargas y contornos para todos los conductores
todas_cargas = []
todos_control = []
Phi = []

for k in range(cantidad_de_conductores):
    cargas_k, control_k = generar_cargas_y_control(cargas_ficticias_por_conductor, radio_conductores, centros_conductores[k],radio_cargas_ficticias)
    todas_cargas.extend(cargas_k)
    todos_control.extend(control_k)
    Phi.extend([potenciales[k]] * cargas_ficticias_por_conductor)
todas_cargas = np.array(todas_cargas)
todos_control = np.array(todos_control)
Phi = np.array(Phi)

#Construir el sistema de ecuaciones del metodo de simulacion de cargas y resolverlo
P = construir_matriz_P(todas_cargas, todos_control)
Q = np.linalg.solve(P, Phi)

#Magnitud de todas las cargas ficticias usadas
print("Magnitudes de todas las cargas ficticias:")
print(Q)
# -----------------------
#METODO DE SIMULACION DE CARGAS
# -----------------------





# -----------------------
#GRAFICAR RESULTADO
# -----------------------
#Define los límites de los ejes y la cantidad de puntos para calcular el potencial escalar electrico en cada uno
x_min = -4  #Minimo eje x
x_max = 4   #Maximo eje x
y_min = -4  #Minimo eje y
y_max = 4   #Maximo eje y
puntos_en_x = 500
puntos_en_y = 500

#Crear malla para graficar el potencial escalar electrico
x = np.linspace(x_min, x_max, puntos_en_x)
y = np.linspace(y_min, y_max, puntos_en_y)
X, Y = np.meshgrid(x, y)
V = calcular_potencial_malla(todas_cargas, Q, X, Y)

#Configuracion escala de colores de la grafica
fig, ax = plt.subplots(figsize=(8, 6))
contour = ax.contourf(X, Y, V, 50, cmap='jet', norm=PowerNorm(gamma=0.5, vmin=np.min(V), vmax=np.max(V))) #Mantener gamma en 0.5, es la mejor escala de colores si se usan 500x500 puntos en la malla
plt.colorbar(contour, label='Potencial escalar electrico [V]')

#Dibujar cada conductor
for centro in centros_conductores:
    circle = plt.Circle(centro, radio_conductores, color='black', fill=False, linestyle='--')
    ax.add_patch(circle)

ax.set_aspect('equal')
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_title(f'Potencial eléctrico con {cantidad_de_conductores} conductor/es {voltaje_fase}V r={radio_conductores}m d=50cm')
plt.grid(True)

#Cambiar los límites de los ejes de la grafica, para mostrar mas de cerca el conductor
plt.xlim(-0.5, 0.5)   # Eje x
plt.ylim(-0.5, 0.5)   # Eje y
plt.show()
# -----------------------
# GRAFICAR RESULTADO
# -----------------------





# -----------------------
#CALCULAR POTENCIAL PROMEDIO EN LA SUPERFICIE DE LOS CONDUCTORES TOMANDO MUCHOS MAS PUNTOS DE CONTROL
# -----------------------
#Nueva cantidad de puntos de control sobre la superficie de cada conductor
multiplicador_puntos_de_control_superficie = 100
n_puntos_superficie = multiplicador_puntos_de_control_superficie * cargas_ficticias_por_conductor

promedios_conductores = []

print("RESULTADOS POR CONDUCTOR")
for idx, centro in enumerate(centros_conductores):
    xc, yc = centro
    ang = np.linspace(0, 2 * np.pi, n_puntos_superficie, endpoint=False)
    puntos_superficie = np.stack([xc + radio_conductores * np.cos(ang),
                                  yc + radio_conductores * np.sin(ang)], axis=1)
    
    potenciales_superficie = []

    for punto in puntos_superficie:
        Vp = 0
        for q, rj in zip(Q, todas_cargas):
            dist = np.linalg.norm(punto - rj)
            if dist != 0:
                Vp += q / (4 * np.pi * epsilon0 * dist)
        potenciales_superficie.append(Vp)
    
    promedio = np.mean(potenciales_superficie)
    error_absoluto = abs(promedio - voltaje_fase)
    error_porcentual = (error_absoluto / voltaje_fase) * 100

    promedios_conductores.append(promedio)

    print(f"Conductor {idx + 1}:")
    print(f"  Potencial promedio = {promedio} V")
    print(f"  Error absoluto     = {error_absoluto} V")
    print(f"  Error porcentual     = {error_porcentual} %")
# -----------------------
#CALCULAR POTENCIAL PROMEDIO EN LA SUPERFICIE DE LOS CONDUCTORES TOMANDO MUCHOS MAS PUNTOS DE CONTROL
# -----------------------





# -----------------------
#VALOR DEL POTENCIAL ESCALAR ELECTRICO EN UN PUNTO ARBITRARIO
# -----------------------
potencial_coordenada_x=300 #Este valor tiene que existir en la malla de puntos que esta calculada
potencial_coordenada_y=280 #Este valor tiene que existir en la malla de puntos que esta calculada
print(f"Potencial escalar electrico en el punto cuya posicion en la matriz V es { (potencial_coordenada_x, potencial_coordenada_y) }: {V[potencial_coordenada_x, potencial_coordenada_y]} V")
# -----------------------
#VALOR DEL POTENCIAL ESCALAR ELECTRICO EN UN PUNTO ARBITRARIO
# -----------------------





# -----------------------
#CALCULO DE LA MAGNITUD DE LA INTENSIDAD DE CAMPO ELECTRICO EN CADA PUNTO DE LA MALLA
# -----------------------
#Calcular pasos en x e y para aproximar la derivada
dx = x[1] - x[0]
dy = y[1] - y[0]

#Calcular gradiente negativo del potencial
dVy, dVx = np.gradient(V, dy, dx)  # Ojo: el primer eje es y, el segundo es x

#Intensidad de campo eléctrico como gradiente negativo del potencial
Ex = -dVx
Ey = -dVy

#Magnitud de la intensidad de campo eléctrico
E_magnitud = np.sqrt(Ex**2 + Ey**2)
# -----------------------
#CALCULO DE LA MAGNITUD DE LA INTENSIDAD DE CAMPO ELECTRICO EN CADA PUNTO DE LA MALLA
# -----------------------





# -----------------------
#MOSTRAR LA MAGNITUD DE LA INTENSIDAD DE CAMPO ELECTRICO EN UN PUNTO ARBITRARIO
# -----------------------
E_coordenada_x=6 #Este valor tiene que existir en la malla de puntos que esta calculada
E_coordenada_y=120 #Este valor tiene que existir en la malla de puntos que esta calculada
print(f"Magnitud de la intensidad de campo electrico en el punto cuya posicion en la matriz E_magnitud es { (E_coordenada_x, E_coordenada_y) }: {E_magnitud[E_coordenada_x, E_coordenada_y]} V/m")
# -----------------------
#MOSTRAR LA MAGNITUD DE LA INTENSIDAD DE CAMPO ELECTRICO EN UN PUNTO ARBITRARIO
# -----------------------





# -----------------------
#CAMBIO EN LA MAGNITUD DE LA INTENSIDAD DE CAMPO ELECTRICO EN CIERTA TRAYECTORIA
# -----------------------
#Ubicacion del punto de inicio de la trayectoria en la matriz E_magnitud
punto_inicio_trayectoria = 251   # Índice de columna donde empieza la trayectoria 
fila_trayectoria = 249          # Índice de la fila que se desea analizar

#Extraer los datos de las matrices para generar la gráfica
x_trayectoria = x[punto_inicio_trayectoria:]
E_trayectoria = E_magnitud[fila_trayectoria, punto_inicio_trayectoria:]

#Graficar
plt.figure(figsize=(8, 5))
plt.plot(x_trayectoria, E_trayectoria, 'o-', color='red', label='Magnitud de la intensidad de campo electrico')
plt.xlabel('x (m)')
plt.ylabel('Magnitud de la intensidad de campo electrico (V/m)')
plt.title(f'Trayectoria horizontal a lo largo del eje x (Ubicacion {punto_inicio_trayectoria,fila_trayectoria} en la matriz "E_magnitud")')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
# -----------------------
#CAMBIO EN LA MAGNITUD DE LA INTENSIDAD DE CAMPO ELECTRICO EN CIERTA TRAYECTORIA
# -----------------------





# -----------------------
#VALORES PARA PODER COPIAR Y GRAFICAR EN OTRA GRAFICA PARA LA SUPERPOSICION DE RESULTADOS
# -----------------------
print("VALORES PARA COPIAR SI SE QUIERE HACER UNA SUPERPOSICION DE GRAFICAS:")
print(E_trayectoria)
# -----------------------
#VALORES PARA PODER COPIAR Y GRAFICAR EN OTRA GRAFICA PARA LA SUPERPOSICION DE RESULTADOS
# -----------------------