## Modelos sencillos FDTD para simulaciones acústicas en MATLAB

Todos los ejemplos solo se han probado en MATLAB 2019a, por lo que en versiones anteriores puede que algunas líneas fallen.

Iré subiendo más modelos poco a poco.

#### FDTD_2D_Basico

Este modelo simula un recinto en dos dimensiones, la excitación es un pulso Ricker e incluye unas líneas simples para emular un objeto rígido en el recinto.

#### FDTD_QRD_Unidimensional_2D

Este modelo simula un recinto de dos dimensiones con un difusor QRD. El difusor se configura en el script según los parámetros de diseño (frecuencia de diseño, número primo generador, etc).
Se puede elegir entre difusor o panel plano, tambien tiene dos excitaciones: pulso Ricker y onda sinusoidal.

Ejemplo con `N = 7` y `fd = 1000`:

<img src="http://jmrplens.com/GitHub_FDTD/QRDpulse.gif" width="30%"></img> <img src="http://jmrplens.com/GitHub_FDTD/QRDcw.gif" width="30%"></img>
