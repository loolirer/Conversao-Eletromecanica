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



NumPontos = 100

intervalos_theta = np.linspace(-29.9, 0, NumPontos)


Wc_real = []
Wc_ideal = []
corrente_real_array = np.zeros(NumPontos)
corrente_ideal_array = np.zeros(NumPontos)
fluxo_conc_array = np.zeros(NumPontos)
B_ar_array = np.zeros(NumPontos)
B_ferro_array = np.linspace(0, 1.816, NumPontos)        # Cria um array linear para B_ferro_array

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


plt.figure(figsize=(14, 7))

# Para o núcleo real
plt.subplot(1, 2, 1)
freal = interp1d(fluxo_conc_array, corrente_real_array, kind='cubic')
fluxo_conc_array_interp = np.linspace(0, 0.431227, 50)
corrente_interp_real = freal(fluxo_conc_array_interp)
plt.plot(corrente_interp_real, fluxo_conc_array_interp, color='red', label='Curva Real')        # Adicionando rótulo para a linha
plt.xlim([0, 25])

plt.xlabel('Corrente (A)')
plt.ylabel('Fluxo de Concentração (mWb)')
plt.title('Curva de λ na Bobina em Função da Corrente Elétrica Real')
plt.grid(True)


# Para o núcleo ideal
plt.subplot(1, 2, 2)
fideal = interp1d(fluxo_conc_array, corrente_ideal_array, kind='cubic')
fluxo_conc_array_interp = np.linspace(0, 0.431227, 50)
corrente_interp_ideal = fideal(fluxo_conc_array_interp)
plt.plot(corrente_interp_ideal, fluxo_conc_array_interp, color='red', label='Curva Ideal')        # Adicionando rótulo para a linha
plt.xlim([0, 15])

plt.xlabel('Corrente (A)')
plt.ylabel('Fluxo de Concentração (mWb)')
plt.title('Curva de λ na Bobina em Função da Corrente Elétrica Ideal')
plt.grid(True)

plt.tight_layout()        # Ajusta automaticamente a disposição dos subplots
plt.legend(loc="best", prop={'size': 10})       # Adicionando legenda para o plot geral
plt.show()


pot_acionamento_ideal = np.max(fluxo_conc_array) * np.max(corrente_ideal_array)
pot_acionamento_real = np.max(fluxo_conc_array) * np.max(corrente_real_array)


Area_energia = sp.integrate.trapezoid(I_intervalo,fluxo_conc_real)        # Energia

# Para o caso ideal há uma aproximacao de 50%, ou seja, a potência mecânica é metade da potencia do controlador.
pot_mecanica_ideal = pot_acionamento_ideal/2
pot_mecanica_real = pot_acionamento_real - Area_energia


razao_ideal = pot_mecanica_ideal/pot_acionamento_ideal
razao_real = (pot_mecanica_real)/(pot_acionamento_real)


print("\n""A razão entre a potência mecânica pelo motor e a potência necessária no caso ideal é: "+ f"{razao_ideal:.3f}")
print("A razão entre a potência mecânica pelo motor e a potência necessária no caso real é: "+ f"{razao_real:.3f}""\n")