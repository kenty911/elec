# elec
elec for wrfout

This is a toy model reproduce lightning for wrfout with conversion rate among(Q_rain, Q_graupel, Q_ice, etc).

## How to conduct
You shuod conduct WRF simulation before this and modify the path in `main.py`.

This requires `poetry` and `make`.

After that compile `c` files to `so` with
`make all`

If you do not have suitable python environment use install python dependency with
`poetry install`

and run with poetry on background
`nohup poetry run python main.py &`

to see the log
`tail -f nohup.out`
