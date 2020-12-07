# Product task maps share domain but not codomain
domain_manifold(task_map::PTMX) = domain_manifold(base_task_map(task_map).map1)

function codomain_manifold(task_map::PTMX)
    base_map = base_task_map(task_map)
    PM{codomain_manifold(base_map.map1), codomain_manifold(base_map.map2)}
end
function codomain_manifolds(task_map::PTMX)
    base_map = base_task_map(task_map)
    codomain_manifold(base_map.map1), codomain_manifold(base_map.map2)
end

choose_chart_emb(xme, task_map::PTM, CM, CN) =
    choose_chart_emb(task_map_emb(xme, task_map, CM, CN), CN)
choose_chart_chart(xm, task_map::PTM, CM, CN) =
    choose_chart_chart(task_map_chart(xm, task_map, CM, CN), CN)

base_task_map(task_map::PTM) = task_map

function task_map_emb(xme, task_map::PTM, CM, CN)
    CN1, CN2 = product_chart_split(CN)
    xn1e = task_map_emb(xme, task_map.map1, CM, CN1)
    xn2e = task_map_emb(xme, task_map.map2, CM, CN2)
    [xn1e; xn2e]
end

function task_map_emb_chart(xme, task_map::PTM, CM, CN)
    CN1, CN2 = product_chart_split(CN)
    xn1 = task_map_emb_chart(xme, task_map.map1, CM, CN1)
    xn2 = task_map_emb_chart(xme, task_map.map2, CM, CN2)
    [xn1; xn2]
end

function task_map_chart_emb(xm, task_map::PTM, CM, CN)
    CN1, CN2 = product_chart_split(CN)
    xn1e = task_map_chart_emb(xm, task_map.map1, CM, CN1)
    xn2e = task_map_chart_emb(xm, task_map.map2, CM, CN2)
    [xn1e; xn2e]
end

function task_map_chart(xm, task_map::PTM, CM, CN)
    CN1, CN2 = product_chart_split(CN)
    xn1 = task_map_chart(xm, task_map.map1, CM, CN1)
    xn2 = task_map_chart(xm, task_map.map2, CM, CN2)
    [xn1; xn2]
end