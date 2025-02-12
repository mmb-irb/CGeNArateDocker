run.sh compiles the code and runs four test simulations, with toml files in /simulations: 1j5n.toml, circular.toml, mito.toml, yeast.toml

/input contains initial structures.
/output will contain MD Results

/build, /build_rel contain compiled files

check.py can be used to check consistent energy results for 1j5n and mito simulations. 
To create references use: "cp ./output/mitoenergy.csv ./output/mitoenergy.8.csv"
"cp ./output/1j5nenergy.csv ./output/1j5nenergy.8.csv"



For new simulations, you need to a toml file in /simulation, with an associated initial structure in /input. It will run together with the previous tests with run.sh