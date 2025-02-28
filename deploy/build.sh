# build OpenGLDepthRenderer
cd ../subs/OpenGLDepthRenderer/
git submodule update --init --recursive
mkdir build
cd build
cmake ..
make
cd ../../../deploy

# install packages with conda
conda env create -n dubinstracking --file ../environment.yaml
