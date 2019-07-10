import glob
import json

    
tex_snippet_1plot = """ 
    \\begin{frame}{%s}
        \\includegraphics[width=0.5\\textwidth]{{"%s"}.pdf}
    \\end{frame}
    """

tex_snippet_2plot = """ 
    \\begin{frame}{%s}
        \\includegraphics[width=0.5\\textwidth]{{"%s"}.pdf}
        \\includegraphics[width=0.5\\textwidth]{{"%s"}.pdf}
    \\end{frame}
    """

tex_snippet_4plot = """ 
    \\begin{frame}{%s}
        \\includegraphics[width=0.25\\textwidth]{{"%s"}.pdf}
        \\includegraphics[width=0.25\\textwidth]{{"%s"}.pdf}
        \\includegraphics[width=0.25\\textwidth]{{"%s"}.pdf}
        \\includegraphics[width=0.25\\textwidth]{{"%s"}.pdf}
    \\end{frame}
    """

def write(f,globString,title):
    plotList = glob.glob(globString)

    for i in range(0, len(plotList), 2):
        if i == len(plotList) - 1:
            f.write(tex_snippet_1plot % (title, plotList[i].replace(".pdf","")))
            print plotList[i]
        else:
            f.write(tex_snippet_2plot % (title, plotList[i].replace(".pdf",""), plotList[i+1].replace(".pdf","")))
            print plotList[i]
            print plotList[i+1]

def writeLine(f, line):
    f.write(line + "\n")

def writeSlide(f, runMap, fileString, variable, eras, title):
    n = len(eras)
    width = 1.0 / float(n)
    x = 0
    dx = int(16 / n)
    writeLine(f, "\\begin{frame}{%s}" % (title))
    for e in eras:
        e_tex = e.replace("_", " ")
        # example: "../histos_DataMC_2016_27_Jun_2019_3/DataMC_Electron_LowDM_dPhi1_2016"
        # example: \includegraphics[width=0.25\textwidth]{{"../histos_DataMC_2016_27_Jun_2019_3/DataMC_Electron_LowDM_dPhi1_2016"}.pdf}
        d = runMap[e]
        name = "../%s/%s_%s" % (d, fileString, e)
        writeLine(f, "\\begin{textblock*}{4cm}(%dcm,2cm)" % x)
        writeLine(f, "\\begin{figure}")
        writeLine(f, "\\centering")
        writeLine(f, "\\textbf{%s}\par\medskip" % e_tex)
        writeLine(f, "\\includegraphics[width=1.0\\textwidth]{{\"%s\"}.pdf}" % (name))
        writeLine(f, "\\caption{%s}" % variable)
        writeLine(f, "\\end{figure}")
        writeLine(f, "\\end{textblock*}")
        x += dx
    writeLine(f, "\\end{frame}")

def main():    
    path = "../histos_DataMC_2016_27_Jun_2019_3/"
    json_file = "plots.json" 
    eras = ["2016", "2017", "2018_AB", "2018_CD"]
    regions = ["LowDM", "HighDM"]
    particles = ["Electron", "Muon", "Photon"]
    j = open(json_file)
    f = open("stack_snippet.tex",'w')
    runMap = json.load(j)
    name = "DataMC_%s_%s_*_%s.pdf"
    variables = ["nj", "ht", "met", "metphi", "dPhi1", "dPhi2", "dPhi3", "dPhi4"]
    # latex versions of variables
    variables_tex = {
                "nj"        : "$N_{jets}$", 
                "ht"        : "$H_T$", 
                "met"       : "$\cancel{E}_T$", 
                "metphi"    : "$\phi_{MET}$", 
                "dPhi1"     : "$\Delta\phi_{1}$", 
                "dPhi2"     : "$\Delta\phi_{2}$", 
                "dPhi3"     : "$\Delta\phi_{3}$", 
                "dPhi4"     : "$\Delta\phi_{4}$",
            }
    
    # example: "../histos_DataMC_2016_27_Jun_2019_3/DataMC_Electron_LowDM_dPhi1_2016"

    for p in particles:
        for r in regions:
            for v in variables:
                v_tex = variables_tex[v]
                writeSlide(f, runMap, "DataMC_%s_%s_%s" % (p, r, v), v_tex, eras, "%s CR %s: %s" % (p, r, v_tex))
    
    #write(f, path+"*Electron*.pdf",   "Electron CR")
    #write(f, path+"*Muon*.pdf",       "Muon CR")
    #write(f, path+"*Photon*.pdf",     "Photon CR")

    f.close()
    j.close()

if __name__ == '__main__':
    main()


