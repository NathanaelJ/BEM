# Author | Nathanael Jenkins
# Description | Imports MATLAB-generated aerofoil data for generating a HAWT blade

import adsk.core, adsk.fusion, adsk.cam, traceback
import io
import os

def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui  = app.userInterface

        product = app.activeProduct
        design = adsk.fusion.Design.cast(product)
        rootComp = design.rootComponent

        # Open files with dialog box
        dlg = ui.createFolderDialog()
        dlg.title = 'Open CSV Files'
        
        dlg.showDialog()
        filenames = os.listdir(dlg.folder)
        # Removes .DS_Store
        if '.DS_Store' in filenames:
            filenames.remove('.DS_Store')

        # Imports each .CSV
        for FOIL in range(len(filenames)):
            app = adsk.core.Application.get()
            ui  = app.userInterface
            
            # Get all components in the active design.
            product = app.activeProduct
            design = adsk.fusion.Design.cast(product)
            title = 'Import Spline csv'
            if not design:
                ui.messageBox('No active Fusion design', title)
                return
            
            # Format filename
            filename = filenames[FOIL]
            parts = filename.split('_') 
            temp = (parts[0])
            temp = temp.split('/')
            temp = temp[-1]
            dist = float(temp)

            with io.open(dlg.folder+'/'+filenames[FOIL], 'r', encoding='utf-8-sig') as f:
            
                # Generate points
                points = adsk.core.ObjectCollection.create()
                line = f.readline()
                data = []
                while line:
                    pntStrArr = line.split(',')
                    for pntStr in pntStrArr:
                        try:
                            data.append(float(pntStr))
                        except:
                            break
                
                    if len(data) >= 3 :
                        point = adsk.core.Point3D.create(data[0], data[1], data[2])
                        points.add(point)
                    line = f.readline()
                    data.clear()            
            if points.count:
                sketch = rootComp.sketches.add(rootComp.xYConstructionPlane)
                sketch.sketchCurves.sketchFittedSplines.add(points)
            else:
                ui.messageBox('No valid points', title)  
        
    except:
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))

