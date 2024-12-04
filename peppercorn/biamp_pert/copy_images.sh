arr="original leaky leaky_shadow occluded"

rm images/*
for i in $arr; do
    echo "`pwd`/$i/plots/"
    cp `pwd`/$i/plots/*.png images/
done
