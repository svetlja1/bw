-v - verbose output
--bits [4/8/16..20...24] - size of table for private key->pubkey creation. The bigger, the faster, but it takes more time to generate (23-6min, 24-12min) and takes more memory on the card.

--inputAddress addr.txt - file with addreses
--inputPhrase test.txt - file with phrases
OR
cat test.txt | ./brainWords --inputIn - reading from pipe

default mode: BTC, single sha256

for ETH:
--eth - single keccak
--camp2 - keccak*2031
^^^ two options could be used together
--ethsha256 - sha256 for ETH

Use example files (testEth.txt, addrEth.txt) to see the difference.

Launch program without arguments for examples.


--root - specifies prefix for your word
--suffix - number o combimnation to start from

You may test it:
-v --eth --root "Hello "  --suffix 708550400 --inputAddress addrEth.txt
It is using part "Hello " as a beginning and then brute-force all possibilities using alphabet:
0123456789 .,'!-ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz

if you use -v parameter, you will see test number, you may later use it to restore work, for example:
-v --eth --root "Hello "  --suffix 708550400 --inputAddress addrEth.txt
So it will start from the given test.
