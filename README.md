# DubinsVisibilityTracking
A Dubins vehilce equipped with a gimballed camera tracks a moving through an urban environment while maitining line of sight to the moving point. 

![Solutions](/readme_image.png)
## Dependencies
 - Ubuntu 20.04.5 LTS
 - Blender 3.6.x
## Build Instructions
 - Run ```./deploy/install_dep.sh```
 - Run ```./build.sh```
## Usage Instructions
- Start database ```docker-compose up -d```
- Activate python environment ```conda activate dubinstracking```
- Usage with ```python -m dubinstracking --help```
- Create volumes ```python -m dubinstracking volumes --out data/environments/paper_small/ --csv data/environments/paper_slow/path.csv --radius 400 --alt 300 --cutoff 2e6```
- Create orbits ```python -m dubinstracking circles --volumes data/environments/paper_slow/volumes.json --out data/environments/paper_slow/ --tspeed 5.0 --uavspeed 20.0 --uavturn 50.0```
- Run Jupyter notebook to create plots ```notebooks/uav_trajectory.ipynb```
## Contact
Collin Hague, chague@charlotte.edu
Artur Wolek, awolek@charlotte.edu