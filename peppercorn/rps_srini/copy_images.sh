arr="0_0 0_5 1 4 9"

rm images/*
for i in $arr; do
    echo "/Users/rajiv/Desktop/PhD/towardsCatalyticCRNs/peppercorn/rps/$i/plots/"
    cp /Users/rajiv/Desktop/PhD/towardsCatalyticCRNs/peppercorn/rps/$i/plots/*.png images/
done
