arr="original leaky occluded leaky_shadow"

mkdir -p images
rm images/*
for i in $arr; do
    echo "`pwd`/$i/plots/"
    cp `pwd`/$i/plots/*.png images/
done
