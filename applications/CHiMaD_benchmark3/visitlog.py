# Visit 2.13.0 log file
ScriptVersion = "2.13.0"
if ScriptVersion != Version():
    print "This script is for VisIt %s. It may not work with version %s" % (ScriptVersion, Version())
ShowAllWindows()
OpenDatabase("/Users/stephendewitt/Documents/workspace/PRISMS_workspace/phaseField_dev/applications/CHiMaD_benchmark3/solution-000000.vtu", 0)
# The UpdateDBPluginInfo RPC is not supported in the VisIt module so it will not be logged.
AddPlot("Volume", "phi", 1, 1)
DrawPlots()
Query("SpatialExtents")
NodePick(coord=(0, 0, 0))
NodePick(coord=(10, 0, 0))
NodePick(coord=(0, 0, 0))
NodePick(coord=(5, 0, 0))
NodePick(coord=(10, 0, 0))
NodePick(coord=(5, 0, 0))
NodePick(coord=(7.5, 0, 0))
NodePick(coord=(10, 0, 0))
NodePick(coord=(7.5, 0, 0))
NodePick(coord=(8.75, 0, 0))
NodePick(coord=(0, 0, 0))
NodePick(coord=(0, 10, 0))
NodePick(coord=(0, 0, 0))
NodePick(coord=(0, 5, 0))
NodePick(coord=(0, 10, 0))
NodePick(coord=(0, 5, 0))
NodePick(coord=(0, 7.5, 0))
NodePick(coord=(0, 10, 0))
NodePick(coord=(0, 7.5, 0))
NodePick(coord=(0, 8.75, 0))
OpenDatabase("/Users/stephendewitt/Documents/workspace/PRISMS_workspace/phaseField_dev/applications/CHiMaD_benchmark3/solution-060000.vtu", 0)
AddPlot("Volume", "phi", 1, 1)
DrawPlots()
NodePick(coord=(0, 0, 0))
NodePick(coord=(10, 0, 0))
NodePick(coord=(20, 0, 0))
NodePick(coord=(30, 0, 0))
NodePick(coord=(40, 0, 0))
NodePick(coord=(50, 0, 0))
NodePick(coord=(60, 0, 0))
NodePick(coord=(50, 0, 0))
NodePick(coord=(55, 0, 0))
NodePick(coord=(50, 0, 0))
NodePick(coord=(52.5, 0, 0))
NodePick(coord=(55, 0, 0))
NodePick(coord=(52.5, 0, 0))
NodePick(coord=(53.75, 0, 0))
NodePick(coord=(55, 0, 0))
NodePick(coord=(0, 0, 0))
NodePick(coord=(0, 10, 0))
NodePick(coord=(0, 20, 0))
NodePick(coord=(0, 30, 0))
NodePick(coord=(0, 40, 0))
NodePick(coord=(0, 50, 0))
NodePick(coord=(0, 60, 0))
NodePick(coord=(0, 50, 0))
NodePick(coord=(0, 55, 0))
NodePick(coord=(0, 50, 0))
NodePick(coord=(0, 52.5, 0))
NodePick(coord=(0, 55, 0))
NodePick(coord=(0, 52.5, 0))
NodePick(coord=(0, 53.75, 0))
NodePick(coord=(0, 55, 0))
OpenDatabase("/Users/stephendewitt/Documents/workspace/PRISMS_workspace/phaseField_dev/applications/CHiMaD_benchmark3/solution-120000.vtu", 0)
AddPlot("Volume", "phi", 1, 1)
DrawPlots()
NodePick(coord=(0, 0, 0))
NodePick(coord=(10, 0, 0))
NodePick(coord=(20, 0, 0))
NodePick(coord=(30, 0, 0))
NodePick(coord=(40, 0, 0))
NodePick(coord=(50, 0, 0))
NodePick(coord=(60, 0, 0))
NodePick(coord=(70, 0, 0))
NodePick(coord=(80, 0, 0))
NodePick(coord=(90, 0, 0))
NodePick(coord=(80, 0, 0))
NodePick(coord=(85, 0, 0))
NodePick(coord=(90, 0, 0))
NodePick(coord=(85, 0, 0))
NodePick(coord=(87.5, 0, 0))
NodePick(coord=(85, 0, 0))
NodePick(coord=(86.25, 0, 0))
NodePick(coord=(87.5, 0, 0))
NodePick(coord=(0, 0, 0))
NodePick(coord=(0, 10, 0))
NodePick(coord=(0, 20, 0))
NodePick(coord=(0, 30, 0))
NodePick(coord=(0, 40, 0))
NodePick(coord=(0, 50, 0))
NodePick(coord=(0, 60, 0))
NodePick(coord=(0, 70, 0))
NodePick(coord=(0, 80, 0))
NodePick(coord=(0, 90, 0))
NodePick(coord=(0, 80, 0))
NodePick(coord=(0, 85, 0))
NodePick(coord=(0, 90, 0))
NodePick(coord=(0, 85, 0))
NodePick(coord=(0, 87.5, 0))
NodePick(coord=(0, 85, 0))
NodePick(coord=(0, 86.25, 0))
NodePick(coord=(0, 87.5, 0))
OpenDatabase("/Users/stephendewitt/Documents/workspace/PRISMS_workspace/phaseField_dev/applications/CHiMaD_benchmark3/solution-180000.vtu", 0)
AddPlot("Volume", "phi", 1, 1)
DrawPlots()
NodePick(coord=(0, 0, 0))
NodePick(coord=(10, 0, 0))
NodePick(coord=(20, 0, 0))
NodePick(coord=(30, 0, 0))
NodePick(coord=(40, 0, 0))
NodePick(coord=(50, 0, 0))
NodePick(coord=(60, 0, 0))
NodePick(coord=(70, 0, 0))
NodePick(coord=(80, 0, 0))
NodePick(coord=(90, 0, 0))
NodePick(coord=(100, 0, 0))
NodePick(coord=(110, 0, 0))
NodePick(coord=(120, 0, 0))
NodePick(coord=(110, 0, 0))
NodePick(coord=(115, 0, 0))
NodePick(coord=(120, 0, 0))
NodePick(coord=(115, 0, 0))
NodePick(coord=(117.5, 0, 0))
NodePick(coord=(120, 0, 0))
NodePick(coord=(117.5, 0, 0))
NodePick(coord=(118.75, 0, 0))
NodePick(coord=(0, 0, 0))
NodePick(coord=(0, 10, 0))
NodePick(coord=(0, 20, 0))
NodePick(coord=(0, 30, 0))
NodePick(coord=(0, 40, 0))
NodePick(coord=(0, 50, 0))
NodePick(coord=(0, 60, 0))
NodePick(coord=(0, 70, 0))
NodePick(coord=(0, 80, 0))
NodePick(coord=(0, 90, 0))
NodePick(coord=(0, 100, 0))
NodePick(coord=(0, 110, 0))
NodePick(coord=(0, 120, 0))
NodePick(coord=(0, 110, 0))
NodePick(coord=(0, 115, 0))
NodePick(coord=(0, 120, 0))
NodePick(coord=(0, 115, 0))
NodePick(coord=(0, 117.5, 0))
NodePick(coord=(0, 120, 0))
NodePick(coord=(0, 117.5, 0))
NodePick(coord=(0, 118.75, 0))
OpenDatabase("/Users/stephendewitt/Documents/workspace/PRISMS_workspace/phaseField_dev/applications/CHiMaD_benchmark3/solution-240000.vtu", 0)
AddPlot("Volume", "phi", 1, 1)
DrawPlots()
NodePick(coord=(0, 0, 0))
NodePick(coord=(10, 0, 0))
NodePick(coord=(20, 0, 0))
NodePick(coord=(30, 0, 0))
NodePick(coord=(40, 0, 0))
NodePick(coord=(50, 0, 0))
NodePick(coord=(60, 0, 0))
NodePick(coord=(70, 0, 0))
NodePick(coord=(80, 0, 0))
NodePick(coord=(90, 0, 0))
NodePick(coord=(100, 0, 0))
NodePick(coord=(110, 0, 0))
NodePick(coord=(120, 0, 0))
NodePick(coord=(130, 0, 0))
NodePick(coord=(140, 0, 0))
NodePick(coord=(150, 0, 0))
NodePick(coord=(140, 0, 0))
NodePick(coord=(145, 0, 0))
NodePick(coord=(150, 0, 0))
NodePick(coord=(145, 0, 0))
NodePick(coord=(147.5, 0, 0))
NodePick(coord=(145, 0, 0))
NodePick(coord=(146.25, 0, 0))
NodePick(coord=(0, 0, 0))
NodePick(coord=(0, 10, 0))
NodePick(coord=(0, 20, 0))
NodePick(coord=(0, 30, 0))
NodePick(coord=(0, 40, 0))
NodePick(coord=(0, 50, 0))
NodePick(coord=(0, 60, 0))
NodePick(coord=(0, 70, 0))
NodePick(coord=(0, 80, 0))
NodePick(coord=(0, 90, 0))
NodePick(coord=(0, 100, 0))
NodePick(coord=(0, 110, 0))
NodePick(coord=(0, 120, 0))
NodePick(coord=(0, 130, 0))
NodePick(coord=(0, 140, 0))
NodePick(coord=(0, 150, 0))
NodePick(coord=(0, 140, 0))
NodePick(coord=(0, 145, 0))
NodePick(coord=(0, 150, 0))
NodePick(coord=(0, 145, 0))
NodePick(coord=(0, 147.5, 0))
NodePick(coord=(0, 145, 0))
NodePick(coord=(0, 146.25, 0))
OpenDatabase("/Users/stephendewitt/Documents/workspace/PRISMS_workspace/phaseField_dev/applications/CHiMaD_benchmark3/solution-300000.vtu", 0)
AddPlot("Volume", "phi", 1, 1)
DrawPlots()
NodePick(coord=(0, 0, 0))
NodePick(coord=(10, 0, 0))
NodePick(coord=(20, 0, 0))
NodePick(coord=(30, 0, 0))
NodePick(coord=(40, 0, 0))
NodePick(coord=(50, 0, 0))
NodePick(coord=(60, 0, 0))
NodePick(coord=(70, 0, 0))
NodePick(coord=(80, 0, 0))
NodePick(coord=(90, 0, 0))
NodePick(coord=(100, 0, 0))
NodePick(coord=(110, 0, 0))
NodePick(coord=(120, 0, 0))
NodePick(coord=(130, 0, 0))
NodePick(coord=(140, 0, 0))
NodePick(coord=(150, 0, 0))
NodePick(coord=(160, 0, 0))
NodePick(coord=(170, 0, 0))
NodePick(coord=(180, 0, 0))
NodePick(coord=(170, 0, 0))
NodePick(coord=(175, 0, 0))
NodePick(coord=(170, 0, 0))
NodePick(coord=(172.5, 0, 0))
NodePick(coord=(175, 0, 0))
NodePick(coord=(172.5, 0, 0))
NodePick(coord=(173.75, 0, 0))
NodePick(coord=(0, 0, 0))
NodePick(coord=(0, 10, 0))
NodePick(coord=(0, 20, 0))
NodePick(coord=(0, 30, 0))
NodePick(coord=(0, 40, 0))
NodePick(coord=(0, 50, 0))
NodePick(coord=(0, 60, 0))
NodePick(coord=(0, 70, 0))
NodePick(coord=(0, 80, 0))
NodePick(coord=(0, 90, 0))
NodePick(coord=(0, 100, 0))
NodePick(coord=(0, 110, 0))
NodePick(coord=(0, 120, 0))
NodePick(coord=(0, 130, 0))
NodePick(coord=(0, 140, 0))
NodePick(coord=(0, 150, 0))
NodePick(coord=(0, 160, 0))
NodePick(coord=(0, 170, 0))
NodePick(coord=(0, 180, 0))
NodePick(coord=(0, 170, 0))
NodePick(coord=(0, 175, 0))
NodePick(coord=(0, 170, 0))
NodePick(coord=(0, 172.5, 0))
NodePick(coord=(0, 175, 0))
NodePick(coord=(0, 172.5, 0))
NodePick(coord=(0, 173.75, 0))
OpenDatabase("/Users/stephendewitt/Documents/workspace/PRISMS_workspace/phaseField_dev/applications/CHiMaD_benchmark3/solution-360000.vtu", 0)
AddPlot("Volume", "phi", 1, 1)
DrawPlots()
NodePick(coord=(0, 0, 0))
NodePick(coord=(10, 0, 0))
NodePick(coord=(20, 0, 0))
NodePick(coord=(30, 0, 0))
NodePick(coord=(40, 0, 0))
NodePick(coord=(50, 0, 0))
NodePick(coord=(60, 0, 0))
NodePick(coord=(70, 0, 0))
NodePick(coord=(80, 0, 0))
NodePick(coord=(90, 0, 0))
NodePick(coord=(100, 0, 0))
NodePick(coord=(110, 0, 0))
NodePick(coord=(120, 0, 0))
NodePick(coord=(130, 0, 0))
NodePick(coord=(140, 0, 0))
NodePick(coord=(150, 0, 0))
NodePick(coord=(160, 0, 0))
NodePick(coord=(170, 0, 0))
NodePick(coord=(180, 0, 0))
NodePick(coord=(190, 0, 0))
NodePick(coord=(200, 0, 0))
NodePick(coord=(190, 0, 0))
NodePick(coord=(195, 0, 0))
NodePick(coord=(200, 0, 0))
NodePick(coord=(195, 0, 0))
NodePick(coord=(197.5, 0, 0))
NodePick(coord=(200, 0, 0))
NodePick(coord=(197.5, 0, 0))
NodePick(coord=(198.75, 0, 0))
NodePick(coord=(0, 0, 0))
NodePick(coord=(0, 10, 0))
NodePick(coord=(0, 20, 0))
NodePick(coord=(0, 30, 0))
NodePick(coord=(0, 40, 0))
NodePick(coord=(0, 50, 0))
NodePick(coord=(0, 60, 0))
NodePick(coord=(0, 70, 0))
NodePick(coord=(0, 80, 0))
NodePick(coord=(0, 90, 0))
NodePick(coord=(0, 100, 0))
NodePick(coord=(0, 110, 0))
NodePick(coord=(0, 120, 0))
NodePick(coord=(0, 130, 0))
NodePick(coord=(0, 140, 0))
NodePick(coord=(0, 150, 0))
NodePick(coord=(0, 160, 0))
NodePick(coord=(0, 170, 0))
NodePick(coord=(0, 180, 0))
NodePick(coord=(0, 190, 0))
NodePick(coord=(0, 200, 0))
NodePick(coord=(0, 190, 0))
NodePick(coord=(0, 195, 0))
NodePick(coord=(0, 200, 0))
NodePick(coord=(0, 195, 0))
NodePick(coord=(0, 197.5, 0))
NodePick(coord=(0, 200, 0))
NodePick(coord=(0, 197.5, 0))
NodePick(coord=(0, 198.75, 0))
OpenDatabase("/Users/stephendewitt/Documents/workspace/PRISMS_workspace/phaseField_dev/applications/CHiMaD_benchmark3/solution-420000.vtu", 0)
AddPlot("Volume", "phi", 1, 1)
DrawPlots()
NodePick(coord=(0, 0, 0))
NodePick(coord=(10, 0, 0))
NodePick(coord=(20, 0, 0))
NodePick(coord=(30, 0, 0))
NodePick(coord=(40, 0, 0))
NodePick(coord=(50, 0, 0))
NodePick(coord=(60, 0, 0))
NodePick(coord=(70, 0, 0))
NodePick(coord=(80, 0, 0))
NodePick(coord=(90, 0, 0))
NodePick(coord=(100, 0, 0))
NodePick(coord=(110, 0, 0))
NodePick(coord=(120, 0, 0))
NodePick(coord=(130, 0, 0))
NodePick(coord=(140, 0, 0))
NodePick(coord=(150, 0, 0))
NodePick(coord=(160, 0, 0))
NodePick(coord=(170, 0, 0))
NodePick(coord=(180, 0, 0))
NodePick(coord=(190, 0, 0))
NodePick(coord=(200, 0, 0))
NodePick(coord=(210, 0, 0))
NodePick(coord=(220, 0, 0))
NodePick(coord=(230, 0, 0))
NodePick(coord=(220, 0, 0))
NodePick(coord=(225, 0, 0))
NodePick(coord=(220, 0, 0))
NodePick(coord=(222.5, 0, 0))
NodePick(coord=(225, 0, 0))
NodePick(coord=(222.5, 0, 0))
NodePick(coord=(223.75, 0, 0))
NodePick(coord=(225, 0, 0))
NodePick(coord=(0, 0, 0))
NodePick(coord=(0, 10, 0))
NodePick(coord=(0, 20, 0))
NodePick(coord=(0, 30, 0))
NodePick(coord=(0, 40, 0))
NodePick(coord=(0, 50, 0))
NodePick(coord=(0, 60, 0))
NodePick(coord=(0, 70, 0))
NodePick(coord=(0, 80, 0))
NodePick(coord=(0, 90, 0))
NodePick(coord=(0, 100, 0))
NodePick(coord=(0, 110, 0))
NodePick(coord=(0, 120, 0))
NodePick(coord=(0, 130, 0))
NodePick(coord=(0, 140, 0))
NodePick(coord=(0, 150, 0))
NodePick(coord=(0, 160, 0))
NodePick(coord=(0, 170, 0))
NodePick(coord=(0, 180, 0))
NodePick(coord=(0, 190, 0))
NodePick(coord=(0, 200, 0))
NodePick(coord=(0, 210, 0))
NodePick(coord=(0, 220, 0))
NodePick(coord=(0, 230, 0))
NodePick(coord=(0, 220, 0))
NodePick(coord=(0, 225, 0))
NodePick(coord=(0, 220, 0))
NodePick(coord=(0, 222.5, 0))
NodePick(coord=(0, 225, 0))
NodePick(coord=(0, 222.5, 0))
NodePick(coord=(0, 223.75, 0))
NodePick(coord=(0, 225, 0))
OpenDatabase("/Users/stephendewitt/Documents/workspace/PRISMS_workspace/phaseField_dev/applications/CHiMaD_benchmark3/solution-480000.vtu", 0)
AddPlot("Volume", "phi", 1, 1)
DrawPlots()
NodePick(coord=(0, 0, 0))
NodePick(coord=(10, 0, 0))
NodePick(coord=(20, 0, 0))
NodePick(coord=(30, 0, 0))
NodePick(coord=(40, 0, 0))
NodePick(coord=(50, 0, 0))
NodePick(coord=(60, 0, 0))
NodePick(coord=(70, 0, 0))
NodePick(coord=(80, 0, 0))
NodePick(coord=(90, 0, 0))
NodePick(coord=(100, 0, 0))
NodePick(coord=(110, 0, 0))
NodePick(coord=(120, 0, 0))
NodePick(coord=(130, 0, 0))
NodePick(coord=(140, 0, 0))
NodePick(coord=(150, 0, 0))
NodePick(coord=(160, 0, 0))
NodePick(coord=(170, 0, 0))
NodePick(coord=(180, 0, 0))
NodePick(coord=(190, 0, 0))
NodePick(coord=(200, 0, 0))
NodePick(coord=(210, 0, 0))
NodePick(coord=(220, 0, 0))
NodePick(coord=(230, 0, 0))
NodePick(coord=(240, 0, 0))
NodePick(coord=(250, 0, 0))
NodePick(coord=(240, 0, 0))
NodePick(coord=(245, 0, 0))
NodePick(coord=(250, 0, 0))
NodePick(coord=(245, 0, 0))
NodePick(coord=(247.5, 0, 0))
NodePick(coord=(250, 0, 0))
NodePick(coord=(247.5, 0, 0))
NodePick(coord=(248.75, 0, 0))
NodePick(coord=(250, 0, 0))
NodePick(coord=(0, 0, 0))
NodePick(coord=(0, 10, 0))
NodePick(coord=(0, 20, 0))
NodePick(coord=(0, 30, 0))
NodePick(coord=(0, 40, 0))
NodePick(coord=(0, 50, 0))
NodePick(coord=(0, 60, 0))
NodePick(coord=(0, 70, 0))
NodePick(coord=(0, 80, 0))
NodePick(coord=(0, 90, 0))
NodePick(coord=(0, 100, 0))
NodePick(coord=(0, 110, 0))
NodePick(coord=(0, 120, 0))
NodePick(coord=(0, 130, 0))
NodePick(coord=(0, 140, 0))
NodePick(coord=(0, 150, 0))
NodePick(coord=(0, 160, 0))
NodePick(coord=(0, 170, 0))
NodePick(coord=(0, 180, 0))
NodePick(coord=(0, 190, 0))
NodePick(coord=(0, 200, 0))
NodePick(coord=(0, 210, 0))
NodePick(coord=(0, 220, 0))
NodePick(coord=(0, 230, 0))
NodePick(coord=(0, 240, 0))
NodePick(coord=(0, 250, 0))
NodePick(coord=(0, 240, 0))
NodePick(coord=(0, 245, 0))
NodePick(coord=(0, 250, 0))
NodePick(coord=(0, 245, 0))
NodePick(coord=(0, 247.5, 0))
NodePick(coord=(0, 250, 0))
NodePick(coord=(0, 247.5, 0))
NodePick(coord=(0, 248.75, 0))
NodePick(coord=(0, 250, 0))
OpenDatabase("/Users/stephendewitt/Documents/workspace/PRISMS_workspace/phaseField_dev/applications/CHiMaD_benchmark3/solution-540000.vtu", 0)
AddPlot("Volume", "phi", 1, 1)
DrawPlots()
NodePick(coord=(0, 0, 0))
NodePick(coord=(10, 0, 0))
NodePick(coord=(20, 0, 0))
NodePick(coord=(30, 0, 0))
NodePick(coord=(40, 0, 0))
NodePick(coord=(50, 0, 0))
NodePick(coord=(60, 0, 0))
NodePick(coord=(70, 0, 0))
NodePick(coord=(80, 0, 0))
NodePick(coord=(90, 0, 0))
NodePick(coord=(100, 0, 0))
NodePick(coord=(110, 0, 0))
NodePick(coord=(120, 0, 0))
NodePick(coord=(130, 0, 0))
NodePick(coord=(140, 0, 0))
NodePick(coord=(150, 0, 0))
NodePick(coord=(160, 0, 0))
NodePick(coord=(170, 0, 0))
NodePick(coord=(180, 0, 0))
NodePick(coord=(190, 0, 0))
NodePick(coord=(200, 0, 0))
NodePick(coord=(210, 0, 0))
NodePick(coord=(220, 0, 0))
NodePick(coord=(230, 0, 0))
NodePick(coord=(240, 0, 0))
NodePick(coord=(250, 0, 0))
NodePick(coord=(260, 0, 0))
NodePick(coord=(270, 0, 0))
NodePick(coord=(280, 0, 0))
NodePick(coord=(270, 0, 0))
NodePick(coord=(275, 0, 0))
NodePick(coord=(270, 0, 0))
NodePick(coord=(272.5, 0, 0))
NodePick(coord=(275, 0, 0))
NodePick(coord=(272.5, 0, 0))
NodePick(coord=(273.75, 0, 0))
NodePick(coord=(0, 0, 0))
NodePick(coord=(0, 10, 0))
NodePick(coord=(0, 20, 0))
NodePick(coord=(0, 30, 0))
NodePick(coord=(0, 40, 0))
NodePick(coord=(0, 50, 0))
NodePick(coord=(0, 60, 0))
NodePick(coord=(0, 70, 0))
NodePick(coord=(0, 80, 0))
NodePick(coord=(0, 90, 0))
NodePick(coord=(0, 100, 0))
NodePick(coord=(0, 110, 0))
NodePick(coord=(0, 120, 0))
NodePick(coord=(0, 130, 0))
NodePick(coord=(0, 140, 0))
NodePick(coord=(0, 150, 0))
NodePick(coord=(0, 160, 0))
NodePick(coord=(0, 170, 0))
NodePick(coord=(0, 180, 0))
NodePick(coord=(0, 190, 0))
NodePick(coord=(0, 200, 0))
NodePick(coord=(0, 210, 0))
NodePick(coord=(0, 220, 0))
NodePick(coord=(0, 230, 0))
NodePick(coord=(0, 240, 0))
NodePick(coord=(0, 250, 0))
NodePick(coord=(0, 260, 0))
NodePick(coord=(0, 270, 0))
NodePick(coord=(0, 280, 0))
NodePick(coord=(0, 270, 0))
NodePick(coord=(0, 275, 0))
NodePick(coord=(0, 270, 0))
NodePick(coord=(0, 272.5, 0))
NodePick(coord=(0, 275, 0))
NodePick(coord=(0, 272.5, 0))
NodePick(coord=(0, 273.75, 0))
OpenDatabase("/Users/stephendewitt/Documents/workspace/PRISMS_workspace/phaseField_dev/applications/CHiMaD_benchmark3/solution-600000.vtu", 0)
AddPlot("Volume", "phi", 1, 1)
DrawPlots()
NodePick(coord=(0, 0, 0))
NodePick(coord=(10, 0, 0))
NodePick(coord=(20, 0, 0))
NodePick(coord=(30, 0, 0))
NodePick(coord=(40, 0, 0))
NodePick(coord=(50, 0, 0))
NodePick(coord=(60, 0, 0))
NodePick(coord=(70, 0, 0))
NodePick(coord=(80, 0, 0))
NodePick(coord=(90, 0, 0))
NodePick(coord=(100, 0, 0))
NodePick(coord=(110, 0, 0))
NodePick(coord=(120, 0, 0))
NodePick(coord=(130, 0, 0))
NodePick(coord=(140, 0, 0))
NodePick(coord=(150, 0, 0))
NodePick(coord=(160, 0, 0))
NodePick(coord=(170, 0, 0))
NodePick(coord=(180, 0, 0))
NodePick(coord=(190, 0, 0))
NodePick(coord=(200, 0, 0))
NodePick(coord=(210, 0, 0))
NodePick(coord=(220, 0, 0))
NodePick(coord=(230, 0, 0))
NodePick(coord=(240, 0, 0))
NodePick(coord=(250, 0, 0))
NodePick(coord=(260, 0, 0))
NodePick(coord=(270, 0, 0))
NodePick(coord=(280, 0, 0))
NodePick(coord=(290, 0, 0))
NodePick(coord=(300, 0, 0))
NodePick(coord=(290, 0, 0))
NodePick(coord=(295, 0, 0))
NodePick(coord=(300, 0, 0))
NodePick(coord=(295, 0, 0))
NodePick(coord=(297.5, 0, 0))
NodePick(coord=(295, 0, 0))
NodePick(coord=(296.25, 0, 0))
NodePick(coord=(0, 0, 0))
NodePick(coord=(0, 10, 0))
NodePick(coord=(0, 20, 0))
NodePick(coord=(0, 30, 0))
NodePick(coord=(0, 40, 0))
NodePick(coord=(0, 50, 0))
NodePick(coord=(0, 60, 0))
NodePick(coord=(0, 70, 0))
NodePick(coord=(0, 80, 0))
NodePick(coord=(0, 90, 0))
NodePick(coord=(0, 100, 0))
NodePick(coord=(0, 110, 0))
NodePick(coord=(0, 120, 0))
NodePick(coord=(0, 130, 0))
NodePick(coord=(0, 140, 0))
NodePick(coord=(0, 150, 0))
NodePick(coord=(0, 160, 0))
NodePick(coord=(0, 170, 0))
NodePick(coord=(0, 180, 0))
NodePick(coord=(0, 190, 0))
NodePick(coord=(0, 200, 0))
NodePick(coord=(0, 210, 0))
NodePick(coord=(0, 220, 0))
NodePick(coord=(0, 230, 0))
NodePick(coord=(0, 240, 0))
NodePick(coord=(0, 250, 0))
NodePick(coord=(0, 260, 0))
NodePick(coord=(0, 270, 0))
NodePick(coord=(0, 280, 0))
NodePick(coord=(0, 290, 0))
NodePick(coord=(0, 300, 0))
NodePick(coord=(0, 290, 0))
NodePick(coord=(0, 295, 0))
NodePick(coord=(0, 300, 0))
NodePick(coord=(0, 295, 0))
NodePick(coord=(0, 297.5, 0))
NodePick(coord=(0, 295, 0))
NodePick(coord=(0, 296.25, 0))
