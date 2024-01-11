arr="0_0 0_5 1 4 9 original leaky occluded leaky_shadow"

rm images/*
for i in $arr; do
    echo "`pwd`/$i/plots/"
    cp `pwd`/$i/plots/*.png images/
done
