#Author-Nathanael Jenkins
#Description-Imports MATLAB-generated aerofoil section data, and lofts to form a bRep blade

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

        # Create loft feature input
        loftFeats = rootComp.features.loftFeatures
        loftInput = loftFeats.createInput(adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        loftSectionsObj = loftInput.loftSections

        dlg = ui.createFolderDialog()
        dlg.title = 'Open CSV Files'
        #if dlg.showOpen() != adsk.core.DialogResults.DialogOK :
        #    return
        
        dlg.showDialog()
        filenames = os.listdir(dlg.folder)
        if '.DS_Store' in filenames:
            filenames.remove('.DS_Store')
        #filenames.sort()

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
            
            filename = filenames[FOIL]
            parts = filename.split('_') 
            temp = (parts[0])
            temp = temp.split('/')
            temp = temp[-1]
            dist = float(temp)

            # Add construction plane with offset
            # ctorPlanes = rootComp.constructionPlanes
            # ctorPlaneInput = ctorPlanes.createInput()
            # offset = adsk.core.ValueInput.createByReal(dist)
            # ctorPlaneInput.setByOffset(rootComp.xYConstructionPlane, offset)
            # ctorPlane1 = ctorPlanes.add(ctorPlaneInput)

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
            profile1 = sketch.profiles.item(0)
            loftSectionsObj.add(profile1)

        # LOFT
        loftInput.isSolid = True
        loftFeats.add(loftInput)
        


    except:
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))

