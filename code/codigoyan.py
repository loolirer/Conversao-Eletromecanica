import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d
from scipy.misc import derivative
from scipy import integrate

import warnings
warnings.filterwarnings('ignore')

# Pontos fornecidos na tabela da página
H = np.array([0, 68, 135, 203, 271, 338, 406, 474, 542, 609, 1100, 1500, 2500, 4000, 5000, 9000, 12000, 20000, 25000])  # Valores de H
B = np.array([0, 0.733, 1.205, 1.424, 1.517, 1.56, 1.588, 1.617, 1.631, 1.646, 1.689, 1.703, 1.724, 1.731, 1.738, 1.761, 1.770, 1.8, 1.816])  # Valores de B

# Interpolação spline
f = interp1d( B, H, kind='cubic')

# Valores de H para os quais queremos interpolar B
B_interp = np.linspace(0, 1.816, 100)
H_interp = f(B_interp)

# Plotando o gráfico BxH

plt.figure(figsize=(14, 7))
plt.plot(H, B, 'o', label='Pontos fornecidos')
plt.plot(H_interp, B_interp, '-', label='Interpolação spline')
plt.xlabel('H (A/m)')
plt.ylabel('B (T)')
plt.title('Curva BxH do Material')
plt.legend()
plt.grid(True)
plt.show()




#Apresente um gráfico do fluxo concatenado na bobina 1 em função da corrente
# aplicada nessa bobina considerando a posição do rotor variando da posição -30°
# até +30 graus.

r = 6.3e-2
D = 8e-2
n = 90
lf = 75e-2 #Professor definiu
g = 2*0.45e-3
u0 = 4*np.pi*1e-7


indice = 0
theta = (np.radians(np.array([-29.9, -20, -10, 0, 10, 20, 29.9])))   # ângulos em radianos
modulo_theta = np.array([])
area_ar = np.array([])

#area do ferro constante
area_ferro = D*r*((np.pi)/6)


#3 matrizes, cada uma com 7 linhas e 100 colunas
flux_conc_ferro = np.empty((7, 100)) #Matriz de arrays
B_ar = np.empty((7, 100)) #np.array([])
corrente = np.empty((7, 100)) #np.array([])
correnteid = np.empty((7,100))


# Enquanto o contador for menor que 7, o bloco de código será executado
while indice < 7:
    modulo_theta = np.append(modulo_theta, [np.abs(theta[indice])])
    #area do ar que varia com o theta
    area_ar = np.append(area_ar, (D*r*(((np.pi)/6)-modulo_theta[indice])))
    #print(area_ar)
    add = 0
    for i in np.arange(0, 1.816, 0.01834): # 1.816 = valor maximo de B
      #fluxo concatenado do ferro
      B_ar[indice][add] = (i*area_ferro)/(area_ar[indice])
      flux_conc_ferro[indice][add] =  ((n*area_ar[indice])*(B_ar[indice][add]))
      #print(lf, n, u0, g)
      #print(H_interp[add])
      corrente[indice][add] = (H_interp[add]*lf + B_ar[indice][add]*g/u0)/n
      correnteid[indice][add]    = (B_ar[indice][add]*g/u0)/n

      add += 1

    # Incrementamos o contador em 1 a cada iteração
    indice += 1

print((B_ar))
#print(flux_conc_ferro)
#print(corrente)


f2 = np.array([])
flux_conc_ferro_interp = np.array([])
corrente_interp = np.array([])

for j in range (0,7):
  f2 = interp1d(flux_conc_ferro[j], corrente[j], kind='cubic')
  f3 = interp1d(flux_conc_ferro[j], correnteid[j], kind = 'cubic')

  # Valores de H para os quais queremos interpolar B
  flux_conc_ferro_interp = np.linspace(0, 0.431227, 50)
  corrente_interp = (f2(flux_conc_ferro_interp))
  corrente_interp_ideal = (f3(flux_conc_ferro_interp))




plt.figure(figsize=(14, 7))

# Para cada ângulo, plotamos o gráfico da corrente real
for j in range(7):
    plt.plot(corrente[j], flux_conc_ferro[j], '-', label=f'Ângulo {np.degrees(theta[j]):.1f}° - Corrente Real')
    plt.xlabel('Corrente (A)')
    plt.ylabel('Fluxo de Concentração de Ferro')
    plt.xlim([-0.1,50])
    plt.title('Fluxo concatenado na Bobina 1 em função da Corrente aplicada')
    plt.legend()
    plt.grid(True)
    plt.show()

# Para cada ângulo, plotamos o gráfico da corrente ideal
for j in range(7):
    plt.plot(correnteid[j], flux_conc_ferro[j], '--', label=f'Ângulo {np.degrees(theta[j]):.1f}° - Corrente Ideal')
    plt.xlabel('Corrente (A)')
    plt.ylabel('Fluxo de Concentração de Ferro')
    plt.xlim([-0.1,50])
    plt.title('Fluxo concatenado na Bobina 1 em função da Corrente aplicada')
    plt.legend()
    plt.grid(True)
    plt.show()

#

H_max = f(1.8)

I_max = (H_max*lf + ((1.8*area_ferro)/(area_ferro)/u0)*(g))/n

print(I_max)

NumPontos = 100

intervalos_theta = np.linspace(-29.9, 29.9, NumPontos)

Wc_real = []
Wc_ideal = []
corrente_real_array = np.zeros(NumPontos)
corrente_ideal_array = np.zeros(NumPontos)
fluxo_conc_array = np.zeros(NumPontos)
B_ar_array = np.zeros(NumPontos)
B_ferro_array = np.linspace(0, 1.816, NumPontos)  # Criar um array linear para B_ferro_array

A = lambda theta: (30 - abs(theta)) * np.pi * r * D / 180

for t in range(len(intervalos_theta)):
  for b in range(len(B_ferro_array)):
    B_ar_array[b] = (B_ferro_array[b] * area_ferro) / (A(intervalos_theta[t]))
    fluxo_conc_array[b] = ((n * A(intervalos_theta[t])) * (B_ar_array[b]))
    corrente_real_array[b] = (H_interp[b] * lf + B_ar_array[b] * g / u0) / n
    corrente_ideal_array[b] = (B_ar_array[b] * g / u0) / n


  fluxo_interpolado_real  = CubicSpline(corrente_real_array, fluxo_conc_array, bc_type='natural')
  fluxo_interpolado_ideal = CubicSpline(corrente_ideal_array, fluxo_conc_array , bc_type='natural')

  # Vetor de corrente para a integração
  I_intervalo = np.linspace(0,I_max,NumPontos)

  fluxo_conc_real  = fluxo_interpolado_real(I_intervalo)
  fluxo_conc_ideal = fluxo_interpolado_ideal(I_intervalo)

  # Definindo o vetor energia integrando o fluxo
  Wc_real.append(sp.integrate.trapezoid(fluxo_conc_real, I_intervalo))
  Wc_ideal.append(sp.integrate.trapezoid(fluxo_conc_ideal, I_intervalo))
  # print(f'Energia real: {sp.integrate.trapezoid(Iint, λ_real)}, Energia ideal: {sp.integrate.trapezoid(Iint, λ_ideal)}')

Wc_real = np.array(Wc_real)
Wc_ideal = np.array(Wc_ideal)

print(len(intervalos_theta))
print(len(Wc_real))

# Determinando a curva através Interpolação
coenergia_interpolado_real = CubicSpline(intervalos_theta, Wc_real, bc_type = "natural")
coenergia_interpolado_ideal = CubicSpline(intervalos_theta, Wc_ideal, bc_type = "natural")

# Achando o Torque a partir da derivada da energia
Torque_real = derivative(coenergia_interpolado_real, intervalos_theta, dx=np.max(intervalos_theta)/NumPontos)
Torque_ideal = derivative(coenergia_interpolado_ideal, intervalos_theta, dx=np.max(intervalos_theta)/NumPontos)


plt.subplot(2,2,1)
plt.plot(intervalos_theta, corrente_real_array)
plt.legend(["Corrente Real"])
plt.subplot(2,2,2)
# plt.plot(theta, I_ideal)
# plt.legend(["Corrente Ideal"])
# plt.subplot(2,2,3)
plt.plot(intervalos_theta, fluxo_conc_array)
plt.legend(["Fluxo Concatenado"])

plt.subplot(2, 2, 2)
plt.plot(intervalos_theta, Torque_real)
plt.title("Torque vs Deslocamento")
plt.xlabel("Teta(Θ°)")
plt.ylabel("Torque, (T)")

plt.subplot(2, 2, 4)
plt.plot(intervalos_theta, Torque_ideal)
plt.xlabel("Teta(Θ°)")
plt.ylabel("Torque, (T)")

plt.subplot(2, 2, 2).set_title("Caso Real - Torque")
plt.subplot(2, 2, 4).set_title("Caso Ideal - Torque")


## Funciona tudo até aqui 


i = 0
incremento = 0.3632
add = 0
flux_conc_ferro = np.array([])
B_ar = np.array([])
corrente = np.array([])
#while com o B do ferro
while i < 1.817:
   #fluxo concatenado do ferro
   flux_conc_ferro = np.append(flux_conc_ferro, (n*area_ferro*i))
   B_ar = np.append(B_ar, ((i*area_ferro)/area_ar[add]))
   corrente = np.append(corrente, (1/n) * [(H_interp*lf) + ((B_ar/u0)*g)]) #pq l ferro esse valor?
   add += 1
   # Incrementamos o contador a cada iteração
   i += incremento


# Interpolação spline
f2 = interp1d(flux_conc_ferro, corrente, kind='cubic')

# Valores de H para os quais queremos interpolar B
flux_conc_ferro_interp = np.linspace(0, 1.816, 50) #Alterar aqui
corrente_interp = f2(flux_conc_ferro_interp)



# Matrizes vazias para armazenar os resultados da interpolação
interpolated_flux_conc_ferro = np.empty_like(flux_conc_ferro)
interpolated_corrente = np.empty_like(corrente)

# Loop para interpolar linha a linha
for i in range(flux_conc_ferro.shape[0]):
    # Interpolação unidimensional para fluxo de concentração de ferro
    f_flux_conc_ferro = interp1d(np.arange(flux_conc_ferro.shape[1]), flux_conc_ferro[i], kind='linear')
    interpolated_flux_conc_ferro[i] = f_flux_conc_ferro(np.linspace(0, flux_conc_ferro.shape[1]-1, 100))

    # Interpolação unidimensional para corrente
    f_corrente = interp1d(np.arange(corrente.shape[1]), corrente[i], kind='linear')
    interpolated_corrente[i] = f_corrente(np.linspace(0, corrente.shape[1]-1, 100))

# Plotagem das curvas interpolaradas
plt.figure(figsize=(10, 6))
for i in range(interpolated_flux_conc_ferro.shape[0]):
    plt.plot(interpolated_flux_conc_ferro[i], interpolated_corrente[i], label=f'Curva {i+1}')

plt.xlabel('Fluxo de Concentração de Ferro')
plt.ylabel('Corrente')
plt.title('Interpolação de Curvas')
plt.legend()
plt.grid(True)
plt.show()