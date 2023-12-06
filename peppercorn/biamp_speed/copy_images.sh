arr="0_0 0_5 1 4 9"

rm images/*
for i in $arr; do
    echo "`pwd`/$i/plots/"
    cp `pwd`/$i/plots/*.png images/
done
