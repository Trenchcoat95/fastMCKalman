
# source $fastMCKalman/fastMCKalman/tests/unitTest.sh


alias helpCat=cat
[[ -x "$(command -v pygmentize)" ]] && alias helpCat="pygmentize -O style=borland,linenos=1 -l bash"
init(){
  cat <<HELP_USAGE | helpCat
  makeData
  makeTest
  analyzeLogs
HELP_USAGE
}

makeData(){
    cat <<EOF >  makeData.sh
#!/bin/bash
    root.exe -n -b -l <<\EOF 2>&1 | tee makeData.log
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/fastMCKalman/aliKalman/test/\"")
    gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
   .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g
    .L $fastMCKalman/fastMCKalman/MC/fastSimulationTest.C+g
    AliPDG::AddParticlesToPdgDataBase();
    AliLog::SetPrintRepetitions(0);
    testTPC(10000,kTRUE);            //setup for the looper development
    .q
EOF
   chmod a+x makeData.sh
   ./makeData.sh
}

makePullTest(){
       cat <<EOF >  makePullTest.sh
#!/bin/bash
    echo fastParticle.root >fastParticle.list
    root.exe -n -b -l <<\EOF 2>&1 | tee makePullTest.log
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/fastMCKalman/aliKalman/test/\"")
    gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
    .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g
    .L $fastMCKalman/fastMCKalman/MC/fastSimulationTest.C+g
     .L $fastMCKalman/fastMCKalman/MC/test_fastSimulation.C+g
     initTreeFast()
     testPulls()
    .q
EOF
   chmod a+x makePullTest.sh
   ./makePullTest.sh
}

analyzeLogs(){
   errors=( "short track" "Too few consecutive points" "Rotation failed" "Propagation failed" "Update failed" "Too big chi2"
   "Correct for material failed" "PropagateToMirrorX failed" )
   errorSources=( "fastParticle::reconstructParticleFull:" "fastParticle::reconstructParticle:" )
   for errorSource in  "${errorSources[@]}"; do
      nErrorstot=$(cat makeData.log | grep -c ${errorSource})
      echo  ${errorSource} ${nErrorstot}
      for error in "${errors[@]}"; do
        nErrors=$(cat makeData.log | grep ${errorSource} | grep -c "${error}")
        echo ${errorSource} ${error} ${nErrors}
      done;
   done
}


init