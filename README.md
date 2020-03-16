## Modelos FDTD sencillos y complejos para simulaciones acústicas en MATLAB

Todos los ejemplos solo se han probado en MATLAB 2019a, por lo que en versiones anteriores puede que algunas líneas fallen.

Iré subiendo más modelos poco a poco.

#### FDTD_2D_Basico

Este modelo simula un recinto en dos dimensiones, la excitación es un pulso Ricker e incluye unas líneas simples para emular un objeto rígido en el recinto.

<img src="http://jmrplens.com/GitHub_FDTD/basico.gif" width="50%"></img> <img src="http://jmrplens.com/GitHub_FDTD/basicorecepcion.jpg" width="43%"></img>

#### FDTD_QRD_Unidimensional_2D

Este modelo simula un recinto de dos dimensiones con un difusor QRD. El difusor se configura en el script según los parámetros de diseño (frecuencia de diseño, número primo generador, etc).
Se puede elegir entre difusor o panel plano, tambien tiene dos excitaciones: pulso Ricker y onda sinusoidal.

Ejemplo con `N = 7` y `fd = 1000`:

<img src="http://jmrplens.com/GitHub_FDTD/QRDpulse.gif" width="30%"></img> <img src="http://jmrplens.com/GitHub_FDTD/QRDcw.gif" width="30%"></img>

#### FDTD_Columna

Este modelo simula una columna de altavoces direccionable. Se puede asignar el número de elementos, la distancia entre ellos, el ángulo de inclinación del haz y dos excitaciones: pulso Ricker y onda sinusoidal.

  <a href="http://jmrplens.com/GitHub_FDTD/ColumnaCOMSOL.zip" target="_blank">_Descargar modelo FEM realizado con COMSOL (con el mismo nivel de personalización que en el script de MATLAB)_</a>

Ejemplos con 32 elementos, 30 grados de inclinación y señal a 2 kHz:

Pulso:

<img src="http://jmrplens.com/GitHub_FDTD/columnapulse.gif" width="50%"></img> <img src="http://jmrplens.com/GitHub_FDTD/colresppulse.jpg" width="43%"></img>

Seno:

<img src="http://jmrplens.com/GitHub_FDTD/columnacw.gif" width="50%"></img> <img src="http://jmrplens.com/GitHub_FDTD/colrespcw.jpg" width="43%"></img>

#### FDTD_Columna_coord

Modelo similar al anterior pero en este caso se focaliza el campo acústico en un punto concreto del espacio que se elija (definido por coordenadas), el resultado es más eficiente en el punto de recepción tal como se puede ver en la señal temporal de los receptores.

<img src="http://jmrplens.com/GitHub_FDTD/columnacwcoor.gif" width="50%"></img> <img src="http://jmrplens.com/GitHub_FDTD/colrespcoor.jpg" width="43%"></img>

Esta focalización se consigue añadiendo un retardo que simula una curvatura de la columa que se agrega al retardo para inclinarla virtualmente:

<img src="http://jmrplens.com/GitHub_FDTD/colarccoor.jpg" width="43%"></img>
