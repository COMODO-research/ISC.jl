function evaluate_Th_SP(Th_SP_value)

    # Make variables global because included scripts use them
    global TotalLength, Gap_BtoB, Edge_Support
    global Th_FL, Th_Web, nt_FL, b_FL, Ht_Sec
    global TW, TG, TD, r, g, searchTol
    global b_SP, Th_SP, nSecHoles
    global Th_WS, Wd_WS, Ht_WS
    global Name1, Name2, Name3, Job_Name
    global BaseDir, Dir_INP, Dir_ODB, Dir_CSV, Dir_PNG

   # INPUTS - IBEAM
    TotalLength         = 70.0                 # Total length               - Beam to Beam
    Gap_BtoB            = 1.0                  # Gap                        - Flange to Flange / Beam to Beam
    Edge_Support        = 50.0                 # Edge Dist. to support      - Beam to Beam
    Th_FL               = 6.0                  # Thickness                  - Flange
    Th_Web              = 6.0                  # Thickness                  - Web
    nt_FL               = 1                    # no. of teeth (each side)   - Flange 
    b_FL                = 50.0                 # Width                      - Flange   
    Ht_Sec              = 50.0                 # Height                     - Section
    TW                  = 10.0                 # Tooth Width                - Flange
    TG                  = 11.0                 # Tooth Gap                  - Flange  
    TD                  = 11.5                 # Tooth Depth                - Flange
    r                   = 1.0                  # Corner Radii               - Flange   
    g                   = 1.0                  # Tolerance                  - b/w Tooth and Holes
    searchTol           = 1e-6                 # Search Tolerance           - General

    # INPUTS - SP
    b_SP                = 20.0                 # Width                      - Side Plate
    Th_SP               = Th_SP_value                  # Thickness                  - Side Plate
    nSecHoles           = 2*nt_FL              # Holes                      - Side Plate

    # INPUTS - WEB STIFFENER
    Th_WS               = Th_Web               # Thickness                  - Web Stiffener 
    Wd_WS               = 0.5*b_FL-0.5*Th_Web  # Width                      - Web Stiffener
    Ht_WS               = Ht_Sec-2*Th_FL       # Height                     - Web Stiffener

    # ABAQUS JOB DETAILS
    ThSP_tag = replace(string(round(Th_SP, digits = 3)), "." => "p")

    Name1    = "TryOpt"
    Name2    = "Th_SP_$ThSP_tag"
    Name3    = "PMcG"
    Job_Name = Name1 * "_" * Name2 * "_" * Name3

    BaseDir             = "C:/Users/23116524/Desktop/Data_Automation"

    Dir_INP             = "$BaseDir/2_INP"
    Dir_ODB             = "$BaseDir/3_ODB"
    Dir_CSV             = "$BaseDir/4_CSV"
    Dir_PNG             = "$BaseDir/5_PNG"

    mkpath(Dir_INP)
    mkpath(Dir_ODB)
    mkpath(Dir_CSV)
    mkpath(Dir_PNG)

    # INCLUDE SCRIPTS
    include("2_1_FL_v1.jl"              )       # B16       ,   A16
    println("Included,        Script, FL, "        , Job_Name)

    include("4_SP_v1.jl"                )       # VSP11     ,   FSP11
    println("Included,        Script, SP, "        , Job_Name)

    # STEP 1 : INP, CREATE
    include("8_1_CreateINP_PMcG_v1.jl"  )
    println("STEP 1,         INP CREATED, "        , Job_Name)

    # STEP 2 : INP,   RUN  
    Abaqus              =  raw"C:\SIMULIA\Commands\abaqus.bat"
    Script_RunINP       =  raw"C:\Users\23116524\Desktop\Automation\1_Scripts\9_RunINP.py"

    println(`$Abaqus cae noGUI=$Script_RunINP -- $Job_Name`  )
    println("Abaqus Job,            Submitted, "   , Job_Name)

    run(`$Abaqus cae noGUI=$Script_RunINP -- $Job_Name`      )
    println("Abaqus Job,            Completed, "   , Job_Name)

    # OPEN ODB, EXPORT CSV
    SMAPython           = raw"C:\SIMULIA\EstProducts\2023\win_b64\code\bin\SMAPython.exe"
    Script_ExportCSV    = raw"C:\Users\23116524\Desktop\Automation\1_Scripts\10_ExportCSV.py"

    println("Abaqus ODB, Open ODB, Export CSV, "   , Job_Name)
    println(`$SMAPython $Script_ExportCSV $Job_Name`         )
    run(`$SMAPython $Script_ExportCSV $Job_Name`             )
    println("Abaqus CSV,           file saved, "   , Job_Name)

    # EXPOER RF2 ~ U2 GRAPH
    include("11_PlotGraph.jl"                                )
    println("Graph saved,     RF2 ~ U2, ", Job_Name          )

    # CONNECT OPTIMISER
    max_rf2 = get_max_rf2(Job_Name, Dir_CSV)

    println("max RF2 : ",max_rf2)

    return max_rf2
end

function get_max_rf2(Job_Name, Dir_CSV)
    csv_path = joinpath(Dir_CSV, Job_Name * ".csv")

    data, header = readdlm(csv_path, ',', header = true)
    header = vec(string.(header))

    rf2_col = findfirst(==("RF2"), header)

    rf2_values = abs.(Float64.(data[:, rf2_col]))

    return maximum(rf2_values)
end