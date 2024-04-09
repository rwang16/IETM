ProjectPath=$(pwd)

Algorithm=IETM
Alpha=1.0
Beta=0.01
# Dataset in {Tweet, SearchSnippets, StackOverflow}
DatasetName=Tweet
# PseudoLongCorpus in {GPT, LLaMa, LLaMa2, DREx}
PseudoLongCorpus=GPT

DataPath=${ProjectPath}/dataset/${DatasetName}
ResultPath=${ProjectPath}/results/${DatasetName}

K=50
for times in {1..5}; do
  java -jar ${ProjectPath}/out/artifacts/IETM_jar/IETM.jar\
       -model ${Algorithm}\
       -dataname ${DatasetName}\
       -alpha ${Alpha}\
       -beta ${Beta}\
       -ntopics ${K}\
       -twords 20\
       -corpus ${DataPath}/${DatasetName}.txt\
       -generateCorpus ${DataPath}/${DatasetName}_${PseudoLongCorpus}.txt\
       -output ${ResultPath}/\
       -name ${Algorithm}_${K}_${Alpha}_${Beta}_${times}\
       -niters 1000
done

# classification evaluation
LabelFile=${DataPath}/${DatasetName}_label.txt
java -jar ${ProjectPath}/out/artifacts/IETM_jar/IETM.jar\
      -model ClassificationEval\
      -label ${LabelFile}\
      -dir ${ResultPath}\
      -prob theta