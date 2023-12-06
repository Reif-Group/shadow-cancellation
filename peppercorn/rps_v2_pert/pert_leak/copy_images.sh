arr="0 1 5 10 20 50 100"
rm images/*
mkdir -p images
for i in $arr; do
    echo "`pwd`/$i/plots/"
    cp `pwd`/$i/plots/*.png images/
done
