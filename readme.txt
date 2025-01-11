For Sale Private SOFT BrainWords $500 $200
Find Lost Bitcoin Passphrases (Brainwallet)
Search passphrases on the fastest program in the world
18,972 passphrases were found, an estimated 5,000 passphrases are lost
They were found on a processor with a maximum speed of 100,000 pasdphrases per second.
Speed RTX 4090 = 360,000,000 passphrases/sec. It's 3600 times faster.
The fact that one 4090 GPU runs in 24 hours is a legendary CPU program brainflayer in 10 years.

To buy the program, telegram https://t.me/cuda8

GPU card	--bits	Speed
4090	24	360 Mkeys/s
A100	24	180 Mkeys/s
A6000	24	180 Mkeys/s
3090	24	180 Mkeys/s
3080 Ti	24	170 Mkeys/s
3080	24	150 Mkeys/s
3070 Ti	24	120 Mkeys/s
3070	24	110 Mkeys/s
3060	24	70 Mkeys/s
2080 S	24	70 Mkeys/s
2070	24	50 Mkeys/s
1 MKey = 1000000 passwords per sec.
How to search for old lost passphrases:
Default alphabet: ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789 .,'!-

-v - verbose output
--bits 16 or 20 or 24
--inputAddress addresses.txt - file with addreses
--inputPhrase dictionary.txt - file with phrases and pass
-d ? (GPU card number, id)
--eth (ETH address search) --inputAddress addrEth.txt default mode: BTC, single sha256

For ETH:
--eth - single keccak
--camp2 - keccak*2031
^^^ two options could be used together
--ethsha256 - sha256 for ETH
--root - specifies prefix for your word
--suffix - number o combimnation to start from

Run: BrainWords.exe -v --eth --root "Hello word " --suffix 708550400 --inputAddress addresses-Eth.txt

Prefix + Combinations
TestA -> TestB -> TestC -> Testzzzzz -> Testzzzzzzzzzzzzzz
Run BrainWords.exe -v --bits 24 --root Test --suffix 0 --inputAddress addresses.txt -d 0

Test test A -> Test test B -> Test test C -> Test test zzzzz -> Test test zzzzzzzzzzzzzz
Run BrainWords.exe -v --bits 24 --root "Test test " --suffix 0 --inputAddress addresses.txt -d 0

Combinations + Suffix (v0.8 only)
A@gmail.com -> B@gmail.com -> C@gmail.com -> zzzzz@gmail.com -> zzzzzzzzzzzzzz@gmail.com
Run BrainWords.exe -v --bits 24 --rootsuffix @gmail.com --suffix 0 --inputAddress addresses.txt -d 0

A Test test -> B Test test -> C Test test -> zzzzz Test test -> zzzzzzzzzzzzzz Test test
Run BrainWords.exe -v --bits 24 --rootsuffix " Test test" --suffix 0 --inputAddress addresses.txt -d 0

Finding passphrases from a text file
Run BrainWords.exe -v --inputPhrase dictionary.txt --inputAddress addresses.txt -d 0

Hashcat + BrainWords (streaming from an external character generator)
Run hashcat.exe --stdout -a 0 dict.txt dict2.txt | BrainWords.exe -v --bits 8 --inputIn --inputAddress addresses.txt -d 0
Or
Run hashcat.exe --stdout -a 0 dict.txt -r use.rule | BrainWords.exe -v --bits 8 --inputIn --inputAddress addresses.txt -d 0
Or
Run ./hashcat -D 2 --stdout -a 3 -i --increment --increment-min=1 --increment-max=8 ?u?l?l?l?d?d?d?d | ./BrainWords -v --bits 8 --inputIn --inputAddress addresses.txt -d 0

Linux:

For ubuntu (linux), be sure to convert the address database to Unix format.
This will remove the ^M from the end of addresses
sudo apt update
sudo apt install -y dos2unix
dos2unix addresses.txt

Run chmod +x BrainWords
Run ./BrainWords -v --bits 24 --root Test --suffix 0 --inputAddress addresses.txt -d 0

Hashcat
Run ./hashcat.bin --stdout -a 6 dictionary.txt ?d?d?d?d | ./BrainWords -v --bits 8 --inputIn --inputAddress addresses.txt -d 0
or
Run ./hashcat.bin --stdout -a 3 --increment ?u?l?l?l?d?d?d | ./BrainWords -v --bits 8 --inputIn --inputAddress addresses.txt -d 0
Low flow rate linux up to 5 Mkeys, Windows up to 10 Mkeys
(If you need more speed, make a copy of the hashcat folder, run )

Frequently asked Questions
How to buy the program?

https://t.me/cuda8

Why did the program freeze at startup?

She didn't hang up! Program start 3090 --bits24 (10 min.)
The program creates tables and downloads to the device
One card requires 4GB or more of RAM to work.
The consumption depends on the size of the table (--bits) and the size of the address file.

What address formats can be uploaded?

BTC bc.., 3.., 1.., or ETH in a text file from a new line
It is recommended to use only OLD addresses 1... from $2

How to continue searching after stopping the program?

The program saves the position
You can start from any position by specifying --suffix 1234567

How to start with 9 characters?

Count the number of characters in the alphabet and choose the one you need.
Example 9 characters (alphabet 64 characters)
64 * 64 * 64 * 64 * 64 * 64 * 64 * 64 = 281474976710656
Use --suffix 281474976710656

How do I search for 9 characters on 50-150 cards?

Divide the desired range into parts (into cards)
18014398509481984 (9) - 281474976710656 (8) = Difference 17732923532771328
17732923532771328 / 50 (cards) = 354658470655426

GPU 0 281474976710656
GPU 1 281474976710656 + 354658470655426 = --suffix 636133447366082
GPU 2 636133447366082 + 354658470655426 = --suffix 990791918021508
GPU 3 990791918021508 + 354658470655426 = --suffix 1345450388676934
...
GPU 50 17378265062115902 + 354658470655426 = --suffix 17732923532771328

Explain what we are looking for? How it works? What's this?

Here is a good example of work for you.
Enter passphrase: fhqyqzhao123 pay attention to the address 1MVFUmYLKmLyC1m3WfyHkEJTZfoHjwDeXE
The difference is that instead of requests to the blockchain.
The program checks against the database of addresses with a positive balance.

What is a brain wallet?

These are passwords or phrases converted to sha256; the output is a private key to the address.
At the address, people stored coins there, a passphrase in their heads.
This method was used 10-15 years ago. At that time, coins were worth almost nothing.
Passwords were lost and forgotten. Do you remember your password 12 years ago? Are you sure?
One character, dot or register is wrong and the entire wallet is lost. There is no restore button.
Look for more information about brain wallet on the Internet. See hire generstion https://brainwalletx.github.io

What does Brain Wallet look like?

Were there any finds at all?

Here is a white list of passwords, phrases, balances, study
https://privatekeyfinder.io/brainwallet/bitcoin/

Multigpu program?

The program is not multigpu. For each card, indicate your id
-d 0 or -d 1 ... -d 11

What does Brain Wallet look like? Were there any finds at all?

Here is a white list of passwords, phrases, balances, study https://privatekeyfinder.io/brainwallet/bitcoin/

Why does the program use a lot of RAM?

The program creates the necessary tables and stores them in memory.

Does the program require an internet connection?

No, the program is looking for the key offline.

The program is sold with source code?

Yes, the program is sold with original source code v0.7 and new v0.8 and instructions

How to change the alphabet in the program?

Yes, open the file Kernel.cuh insert your alphabet.
The program only accepts en + numbers and symbols

#define ALPHABET_LEN 68
const char ALPHABET[69] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789 .,'!-";
device constant char _ALPHABET[69] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789 .,'!-";
ALPHABET_LEN the exact number of characters in the alphabet
_ALPHABET[the exact number of characters in the alphabet + 1]

#define ALPHABET_LEN 36
const char ALPHABET[37] = "0123456789abcdefghijklmnopqrstuvwxyz";
device constant char _ALPHABET[37] = "0123456789abcdefghijklmnopqrstuvwxyz";

How to compile a program?

ubuntu 20.04 CUDA 11.7 (for RTX 4090 CUDA 12.4)
Run: make

vast.ai
Template: cuda:12.0.1-devel-ubuntu20.04

Image: nvidia/cuda:12.0.1-devel-ubuntu20.04

Image CUDA version:
Incompatible images hidden
Launch Type: jupyter

Windows use VS2019 + install CUDA 11.7 (for RTX 4090 CUDA 12.4)

How can I make sure that the program does not stop after it finds it?

Open the main.cu file using a text editor
Line 1100
return true; -> //return true;
Line 1124
uncompResult = true; -> //uncompResult = true;
Line 1140
compResult = true; -> //compResult = true;
Line 1155
bech32Result = true; -> //bech32Result = true;
Line 1170
compResult = true; -> //compResult = true;

In what modes is the program looking for?

She is looking for everyone at once.
Addresses 1 (uncompressed + compressed), addresses 3... in bc...

I have a RTX 3060 TI card, and I have a low speed, how can speed up?
In the new drivers for 30xx Ti, a limiter is installed that slows down the speed by half.
You need to download the old driver from six months ago. 496.13
Delete the new driver, install the old driver, the speed will increase x2
After searching, you can install new drivers.

What arguments are there in the program?

-v Display the generation position in the program window
-b Number of gpu blocks (default set automatically)
-t Number of gpu cores (default set automatically)
--bits Number of bits to generate table 8, 16, 20, 24
--fstatus File name status (default fileStatus.txt)
--inputAddress List of addresses of BTC or ETH in a tex file from a new line
--inputPhrase Name of a text file with a list of BTC or ETH addresses starting on a new line
--inputIn Receive a stream of phrases from an external generator
--iteration Number of iterations SHA256
--eth Enable mode ETH (default BTC)
--camp2 Enable search mode for camp2 (keccak*2031) ETH addresses
--ethsha256 sha256 for ETH
--root Fixed part in generation --root 1234 or --root " 1234 1234" with space
--suffix Number o combimnation to start from --suffix 0
-d number of the required gpu -d 0 or -d 1

