#! /bin/bash
#copies files from the development folder to the standalone folder
cp -R /home/phac/workspace/Panseq_dev/Panseq2/lib/Blast /home/phac/workspace/Panseq_standalone/lib
cp -R /home/phac/workspace/Panseq_dev/Panseq2/lib/FileInteraction /home/phac/workspace/Panseq_standalone/lib
cp -R /home/phac/workspace/Panseq_dev/Panseq2/lib/LociSelector /home/phac/workspace/Panseq_standalone/lib/
cp -R /home/phac/workspace/Panseq_dev/Panseq2/lib/Logging /home/phac/workspace/Panseq_standalone/lib
cp -R /home/phac/workspace/Panseq_dev/Panseq2/lib/MSA/BlastBased /home/phac/workspace/Panseq_standalone/lib/MSA
cp -R /home/phac/workspace/Panseq_dev/Panseq2/lib/Mummer /home/phac/workspace/Panseq_standalone/lib
cp -R /home/phac/workspace/Panseq_dev/Panseq2/lib/Muscle /home/phac/workspace/Panseq_standalone/lib
cp -R /home/phac/workspace/Panseq_dev/Panseq2/lib/NovelRegion /home/phac/workspace/Panseq_standalone/lib
cp -R /home/phac/workspace/Panseq_dev/Panseq2/lib/Pipeline /home/phac/workspace/Panseq_standalone/lib
cp /home/phac/workspace/Panseq_dev/Panseq2/lib/core_config.txt /home/phac/workspace/Panseq_standalone/lib
cp /home/phac/workspace/Panseq_dev/Panseq2/lib/novel_config.txt /home/phac/workspace/Panseq_standalone/lib
cp /home/phac/workspace/Panseq_dev/Panseq2/lib/novelRegionFinder.pl /home/phac/workspace/Panseq_standalone/lib
cp /home/phac/workspace/Panseq_dev/Panseq2/lib/coreAccessory.pl /home/phac/workspace/Panseq_standalone/lib/coreAccessory.pl
echo "Panseq_standalone updated."