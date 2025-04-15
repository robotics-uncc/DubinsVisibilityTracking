# DubinsVisibilityTracking
A Dubins vehilce equipped with a gimballed camera tracks a moving point through an urban environment while maintaining line of sight to the moving point. 

![Solutions](/readme_image.png)
## Dependencies
 - Ubuntu 20.04.5 LTS
 - Blender 3.6.x
 - conda
## Build Instructions
 - Run ```./deploy/install_dep.sh```
 - Run ```cd ./deploy/```
 - Run ```./build.sh```
## Usage Instructions
- Activate python environment ```conda activate dubinstracking```
- Usage with ```python -m dubinstracking --help```
- Create volumes ```python -m dubinstracking volumes --out data/environments/paper_slow/ --csv data/environments/paper_slow/path.csv --radius 400 --alt 300 --cutoff 20e6 --map data/uptownCharlotte.obj```
- Create orbits ```python -m dubinstracking circles --volumes data/environments/paper_slow/volumes.json --out data/environments/paper_slow/ --tspeed 5.0 --uavspeed 20.0 --uavturn 50.0```
- Run Jupyter notebook to create plots ```notebooks/uav_trajectory.ipynb```
## Contact
Collin Hague, chague@charlotte.edu
Artur Wolek, awolek@charlotte.edu
