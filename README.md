# MCOC2020-P1
**Entrega 1**


![Screenshot](https://github.com/Alberto-Hurtado/MCOC2020-P1/blob/master/Entrega%201/balistica.png)


**Entrega 2**

A pesar de los extensos intentos en los calculos, debo tener algún parametro incorrecto o algún error en el planteamiento de las ecuaciones que no logré identificar,
por lo que no pude llegar a una solución para el caso planteado. En los parametros analizados, el satelite siempre llegó a dos destinos: diverger de la orbita o estrellarse en la superficie. Adjunto aqui un aproximado de la nuve de resultados obtenidos.
![Entrega_2_intento_0](https://github.com/Alberto-Hurtado/MCOC2020-P1/blob/master/Entrega%202/Intento_0.png)


**Entrega 5**

1.- 


![orbita_sim_sin_J](https://github.com/Alberto-Hurtado/MCOC2020-P1/blob/master/Entrega%205/orbita_real_h_vs_t_sin_J.png)


2.- Como se puede apreciar en las siguientes imagenes, la deriva del metodo de eulint es considerablemente mayor a la de odeint, teniendo una diferencia de deriva de 20 855 Km. El tiempo que demora cada uno de los metodos es de ~ 0.5 s para eulint y ~ 0.2 s para odeint.


![deriva_odeint](https://github.com/Alberto-Hurtado/MCOC2020-P1/blob/master/Entrega%205/deriva_odeint.png)
![deriva_eulint](https://github.com/Alberto-Hurtado/MCOC2020-P1/blob/master/Entrega%205/deriva_eulint.png)


3.- Después de un Nsubdiviciones = 10 001, se obtuvo una disminucion del error de un ~ 3800 %, pero aun se obtuvo un error del ~ 220 %. Con respecto al tiempo, para Nsubdivisiones=1, el tiempo transcurrido fue de ~ 0.5 s, en cambio para Nsubdivisiones=10 001, el tiempo fue de ~ 1.5 hrs, lo que es un incremento considerable.


![deriva_eulint_N10001](https://github.com/Alberto-Hurtado/MCOC2020-P1/blob/master/Entrega%205/Nsubdivisiones%3D10001%20t%3D4994%2C35%20s.png)


4.- El cambio temporal de hacer los arreglos con J2 y J3 es de ~ 0.1 s, por lo que en teoría aumentó su demanda en un ~ 50 %. Al implementar los terminos de J2 y J3, hubo una mejoría del ~50 % tambien en relacion a la precision. 


![deriva_odeint](https://github.com/Alberto-Hurtado/MCOC2020-P1/blob/master/Entrega%205/deriva_odeint_con_ajuste.png)
