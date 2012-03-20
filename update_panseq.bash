#! /bin/bash
#copies files from the development folder to the standalone folder
cp -R /home/phac/workspace/Panseq_dev/Panseq2/lib/Blast /home/phac/workspace/Panseq_standalone/lib
cp -R /home/phac/workspace/Panseq_dev/Panseq2/lib/FileInteraction /home/phac/workspace/Panseq_standalone/lib
cp -R /home/phac/workspace/Panseq_dev/Panseq2/lib/LociSelector /home/phac/workspace/Panseq_standalone/lib/
cp -R /home/phac/workspace/Panseq_dev/Panseq2/lib/Logging /home/phac/workspace/Panseq_standalone/lib
cp -R /home/phac/workspace/Panseq_dev/Panseq2/lib/MSA/BlastBased /home/phac/workspace/Panseq_standalone/lib
cp -R /home/phac/workspace/Panseq_dev/Panseq2/lib/Mummer /home/phac/workspace/Panseq_standalone/lib
cp -R /home/phac/workspace/Panseq_dev/Panseq2/lib/Muscle /home/phac/workspace/Panseq_standalone/lib
cp -R /home/phac/workspace/Panseq_dev/Panseq2/lib/NovelRegion /home/phac/workspace/Panseq_standalone/lib
cp -R /home/phac/workspace/Panseq_dev/Panseq2/lib/Pipeline /home/phac/workspace/Panseq_standalone/lib
echo "Panseq_standalone updated."