var FILE_ARR_CYCLE = REPLACE_FILE_ARR;
var dark_red_style = 'mapbox://styles/richardsowers/cjohs10mm1f9m2sltnluju4fm';
var street_style = 'mapbox://styles/mapbox/streets-v11';
var basic_style = 'mapbox://styles/richardsowers/cjohy2ivb1jz52so51anbsbq2';
var VZ_style = 'mapbox://styles/mapbox/dark-v9';

function readJSONcycle() {
    Promise.all(FILE_ARR_CYCLE.map(url => fetch(url).then(resp => resp.json())))
        .then(jsons => {
                mapboxgl.accessToken='YOUR API KEY';
               var map = new mapboxgl.Map({
                       container: 'mapEle',
                       style: VZ_style,
                       center: [-73.9725, 40.752],
                       zoom: 12,
                   });
               
               map.on('load', function() {
                       for (const [index,source] of jsons.entries()) {
                           map.addSource('src'+index, {
                               'type':'geojson',
                                   'data': source});

                       map.addLayer({
                               'id': 'src'+index,
                                   'type': 'line',
                                   'source': 'src'+index,
                                   'paint': {
                                   'line-color': ['get','color'],
                                       'line-width': 3
                                       }});
                       }

                       var waitDiv = document.createElement("div");
                       waitDiv.setAttribute("id","waitDivCycle");
                       document.body.appendChild(waitDiv);
                   });
               
            })}
