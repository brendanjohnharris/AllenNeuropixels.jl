
function mouseconnectivitycache()
    mcc = mouse_connectivity_cache.MouseConnectivityCache(manifest_file=mouseconnectivitymanifest)
end
export MouseConnectivityCache


function getannotationvolume()
    vol, info = mouseconnectivitycache().get_annotation_volume()
end

function getstructuremask(structureid::Number)
    strucutreid = Int(structureid)
    mask, info =  mouseconnectivitycache().get_structure_mask(structureid)
end


function gettemplatevolume()
    template, template_info =  mouseconnectivitycache().get_template_volume()
end

function getreferencespace()
    reference, reference_info = mouseconnectivitycache().get_reference_space()
end

