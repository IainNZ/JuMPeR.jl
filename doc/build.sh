# JuMPeR's documentation is built with Judo
# https://github.com/dcjones/Judo.jl
# Install it with:
# cd ~/.julia/v0.3/
# git clone https://github.com/dcjones/Judo.jl.git Judo
# and then keep installing dependencies until it works :D
rm -rf ./html
julia -e 'using Judo; Judo.collate("JuMPeR", template=Pkg.dir("JuMPeR", "doc", "template"))'
cp logo.svg ./html/