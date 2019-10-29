# ,pdf,PDF_Weight,pdfWeight,pdfWeight_Up,pdfWeight_Down,,,,,
# ,pu,PU_Weight,puWeight,puWeight_Up,puWeight_Down,,,,,
# ,isr,ISR_Weight,ISRWeight,ISRWeight_Up,ISRWeight_Down,,,,,
# ,prefire,Prefire_Weight,PrefireWeight,PrefireWeight_Up,PrefireWeight_Down,,,,,
# ,btag,b,BTagWeight,BTagWeight_Up,BTagWeight_Down,,,,,
# ,dielec_sf,dielec_sf,DiElecSF,DiElecSF_Up,DiElecSF_Down,,,,,
# ,dimu_sf,dimu_sf,DiMuSF,DiMuSF_Up,DiMuSF_Down,,,,,
# ,photon_sf,photon_sf,photonSF,photonSF_Up,photonSF_Up,,,,,
# ,dielec_trig,dielec_trig,DiElecTriggerEffPt,DiElecTriggerEffPt_Up,DiElecTriggerEffPt_Down,,,,,
# ,dimu_trig,dimu_trig,DiMuTriggerEffPt,DiMuTriggerEffPt_Up,DiMuTriggerEffPt_Down,,,,,
# ,photon_trig,photon_trig,Stop0l_trigger_eff_Photon_pt,Stop0l_trigger_eff_Photon_pt_up,Stop0l_trigger_eff_Photon_pt_down,,,,,
# ,met_trig_lowdm,trigger_err_lowdm,Stop0l_trigger_eff_MET_low_dm,Stop0l_trigger_eff_MET_low_dm_up,Stop0l_trigger_eff_MET_low_dm_down,,,,,
# ,met_trig_highdm,trigger_err_highdm,Stop0l_trigger_eff_MET_high_dm,Stop0l_trigger_eff_MET_high_dm_up,Stop0l_trigger_eff_MET_high_dm_down,,,,,

# Generate systematics json file from csv

import csv

with open('systematicsJsonGen.py', 'r') as csvfile:
    csv_object = csv.reader(csvfile, delimiter = ',')

    with open('systematics.json', 'w') as jsonFile:
        jsonFile.write("{\n")
        for row in csv_object:
            if (len(row) < 4):
                break
            jsonFile.write('    "' + row[1] + '" : {\n')
            jsonFile.write('          "name"    : "' + row[2] + '",\n')
            jsonFile.write('          "nominal" : "' + row[3] + '",\n')
            jsonFile.write('          "up"      : "' + row[4] + '",\n')
            jsonFile.write('          "down"    : "' + row[5] + '"\n')
            jsonFile.write('    },\n' )
        jsonFile.write( '}\n' )
    
    
    
