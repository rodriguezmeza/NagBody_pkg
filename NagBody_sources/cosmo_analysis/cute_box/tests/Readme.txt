
In data.dat there is a configuration space sample. Can be visualized with

$ nplot2d in=data.dat plotjoined=0  ws=1

Correlation can be done:

$ cute_box params.txt

For a gadget input (check the params_cola_064.txt for details):

$ cute_box params_cola_064.txt
$ nplot2d in=corr128_cola_064.dat