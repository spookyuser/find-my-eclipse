mapboxgl.accessToken = 'pk.eyJ1IjoiY2FsZWJkYXNkZmEiLCJhIjoiY21rY3BhbWsxMDNxNjNnczNtc3BuOHE5cCJ9.bX4Try-h9qXR7p2QpYueqA'

const map = new mapboxgl.Map({
    container: 'map',
    style: 'mapbox://styles/mapbox/satellite-v9',
    projection: 'globe',
    zoom: 1.5,
    center: [0, 20]
});

map.on('style.load', () => {
    map.setFog({
        color: 'rgb(186, 210, 235)',
        'high-color': 'rgb(36, 92, 223)',
        'horizon-blend': 0.02,
        'space-color': 'rgb(11, 11, 25)',
        'star-intensity': 0.6
    });
});

let selectedMarker = null;

function selectLocation(lng, lat, flyTo = true) {
    if (selectedMarker) {
        selectedMarker.remove();
    }

    selectedMarker = new mapboxgl.Marker({ color: '#f5a623' })
        .setLngLat([lng, lat])
        .addTo(map);

    if (flyTo) {
        map.flyTo({ center: [lng, lat], zoom: 8 });
    }

    // Show info panel and hide instructions
    document.getElementById('info-panel').classList.add('visible');
    const instructions = document.getElementById('instructions');
    if (instructions) {
        instructions.style.display = 'none';
    }

    document.getElementById('content').innerHTML = '<p class="loading">Loading...</p>';

    fetch(`/api/eclipses?lat=${lat.toFixed(6)}&lon=${lng.toFixed(6)}`)
        .then(r => r.json())
        .then(data => {
            displayResults(data);
        })
        .catch(err => {
            document.getElementById('content').innerHTML = `<p class="no-eclipse">Error: ${err.message}</p>`;
        });
}

async function searchPlace(query) {
    const response = await fetch(
        `https://api.mapbox.com/geocoding/v5/mapbox.places/${encodeURIComponent(query)}.json?access_token=${mapboxgl.accessToken}&limit=1`
    );
    const data = await response.json();
    if (data.features && data.features.length > 0) {
        const [lng, lat] = data.features[0].center;
        selectLocation(lng, lat, true);
    }
}

// Autocomplete functionality
const searchInput = document.getElementById('search-input');
const autocompleteList = document.getElementById('autocomplete-list');
let autocompleteDebounce = null;
let activeIndex = -1;
let currentSuggestions = [];

async function fetchSuggestions(query) {
    if (!query || query.length < 2) {
        hideSuggestions();
        return;
    }

    try {
        const response = await fetch(
            `https://api.mapbox.com/geocoding/v5/mapbox.places/${encodeURIComponent(query)}.json?access_token=${mapboxgl.accessToken}&limit=5&autocomplete=true`
        );
        const data = await response.json();

        if (data.features && data.features.length > 0) {
            currentSuggestions = data.features;
            showSuggestions(data.features);
        } else {
            hideSuggestions();
        }
    } catch (err) {
        console.error('Autocomplete error:', err);
        hideSuggestions();
    }
}

function showSuggestions(features) {
    autocompleteList.innerHTML = features.map((feature, index) => {
        const name = feature.text || feature.place_name;
        const context = feature.place_name.replace(name + ', ', '').replace(name, '');
        return `
            <div class="autocomplete-item" data-index="${index}">
                <div class="autocomplete-item-name">${escapeHtml(name)}</div>
                ${context ? `<div class="autocomplete-item-context">${escapeHtml(context)}</div>` : ''}
            </div>
        `;
    }).join('');

    autocompleteList.classList.add('visible');
    activeIndex = -1;
}

function hideSuggestions() {
    autocompleteList.classList.remove('visible');
    currentSuggestions = [];
    activeIndex = -1;
}

function escapeHtml(text) {
    const div = document.createElement('div');
    div.textContent = text;
    return div.innerHTML;
}

function selectSuggestion(index) {
    if (currentSuggestions[index]) {
        const feature = currentSuggestions[index];
        const [lng, lat] = feature.center;
        searchInput.value = feature.place_name;
        hideSuggestions();
        selectLocation(lng, lat, true);
    }
}

function updateActiveItem() {
    const items = autocompleteList.querySelectorAll('.autocomplete-item');
    items.forEach((item, i) => {
        item.classList.toggle('active', i === activeIndex);
    });
    if (activeIndex >= 0 && items[activeIndex]) {
        items[activeIndex].scrollIntoView({ block: 'nearest' });
    }
}

searchInput.addEventListener('input', (e) => {
    clearTimeout(autocompleteDebounce);
    autocompleteDebounce = setTimeout(() => {
        fetchSuggestions(e.target.value.trim());
    }, 200);
});

searchInput.addEventListener('keydown', (e) => {
    const items = autocompleteList.querySelectorAll('.autocomplete-item');

    if (e.key === 'ArrowDown') {
        e.preventDefault();
        if (autocompleteList.classList.contains('visible')) {
            activeIndex = Math.min(activeIndex + 1, items.length - 1);
            updateActiveItem();
        }
    } else if (e.key === 'ArrowUp') {
        e.preventDefault();
        if (autocompleteList.classList.contains('visible')) {
            activeIndex = Math.max(activeIndex - 1, 0);
            updateActiveItem();
        }
    } else if (e.key === 'Enter') {
        e.preventDefault();
        if (activeIndex >= 0 && autocompleteList.classList.contains('visible')) {
            selectSuggestion(activeIndex);
        } else {
            const query = searchInput.value.trim();
            if (query) {
                hideSuggestions();
                searchPlace(query);
            }
        }
    } else if (e.key === 'Escape') {
        hideSuggestions();
    }
});

searchInput.addEventListener('blur', () => {
    // Delay to allow click on suggestion
    setTimeout(hideSuggestions, 150);
});

searchInput.addEventListener('focus', () => {
    const query = searchInput.value.trim();
    if (query.length >= 2 && currentSuggestions.length > 0) {
        autocompleteList.classList.add('visible');
    }
});

autocompleteList.addEventListener('click', (e) => {
    const item = e.target.closest('.autocomplete-item');
    if (item) {
        const index = parseInt(item.dataset.index, 10);
        selectSuggestion(index);
        hideSuggestions();
    }
});

document.getElementById('search-btn').addEventListener('click', () => {
    const query = searchInput.value.trim();
    if (query) {
        hideSuggestions();
        searchPlace(query);
    }
});

document.querySelectorAll('.landmark-btn').forEach(btn => {
    btn.addEventListener('click', () => {
        const lng = parseFloat(btn.dataset.lng);
        const lat = parseFloat(btn.dataset.lat);
        selectLocation(lng, lat, true);
    });
});

map.on('click', (e) => {
    const { lng, lat } = e.lngLat;
    selectLocation(lng, lat, false);
});

function formatDuration(seconds) {
    if (!seconds) return null;
    const mins = Math.floor(seconds / 60);
    const secs = Math.round(seconds % 60);
    if (mins > 0) {
        return `${mins}m ${secs}s of totality`;
    }
    return `${secs}s of totality`;
}

function getNasaUrl(dateIso) {
    const ecl = dateIso.replace(/-/g, '');
    return `https://eclipse.gsfc.nasa.gov/SEsearch/SEsearchmap.php?Ecl=${ecl}`;
}

// Layer IDs for our custom eclipse path rendering
const ECLIPSE_LAYER_IDS = [
    'path-fill',
    'path-fill-glow',
    'centerline',
    'centerline-glow',
    'northern-limit',
    'southern-limit'
];

function clearEclipseLayers(prefix) {
    ECLIPSE_LAYER_IDS.forEach(id => {
        const layerId = `${prefix}-${id}`;
        if (map.getLayer(layerId)) {
            map.removeLayer(layerId);
        }
    });
    // Also remove sources
    ['path-polygon', 'centerline', 'northern-limit', 'southern-limit'].forEach(id => {
        const sourceId = `${prefix}-${id}`;
        if (map.getSource(sourceId)) {
            map.removeSource(sourceId);
        }
    });
}

function formatEclipseId(dateIso) {
    // Convert ISO date (2024-04-08) to NASA format (+20240408)
    return '+' + dateIso.replace(/-/g, '');
}

function parseNasaCoords(jsCode, varName) {
    // Extract coordinate arrays from NASA's JS code
    // Format: const varName = [{lat: ..., lng: ...}, {lat: ..., lng: ...}, ...];
    const regex = new RegExp(`(?:const|var)\\s+${varName}\\s*=\\s*\\[([\\s\\S]*?)\\];`);
    const match = jsCode.match(regex);
    if (!match) return null;

    try {
        // Parse the object array format: {lat: 25.5, lng: -105.2}
        const coordsText = match[1];
        const coords = [];
        const pointRegex = /\{lat:\s*([-\d.]+),\s*lng:\s*([-\d.]+)\}/g;
        let pointMatch;
        while ((pointMatch = pointRegex.exec(coordsText)) !== null) {
            const lat = parseFloat(pointMatch[1]);
            const lng = parseFloat(pointMatch[2]);
            coords.push([lng, lat]); // GeoJSON is [lng, lat]
        }
        return coords.length > 0 ? coords : null;
    } catch (e) {
        console.warn(`Failed to parse ${varName}:`, e);
        return null;
    }
}

function parseNasaTimeData(jsCode) {
    // Extract time labels from NASA's code
    // Format: const centralLimitLabels = [{coord: {lat: ..., lng: ...}, label: "time"}, ...];
    const regex = /(?:const|var)\s+centralLimitLabels\s*=\s*\[([\s\S]*?)\];/;
    const match = jsCode.match(regex);
    if (!match) return null;

    try {
        const labelsText = match[1];
        const labels = [];
        // Match: {coord: {lat: 25.5, lng: -105.2}, label: "18:42"}
        const labelRegex = /\{coord:\s*\{lat:\s*([-\d.]+),\s*lng:\s*([-\d.]+)\},\s*label:\s*["']([^"']+)["']\}/g;
        let labelMatch;
        while ((labelMatch = labelRegex.exec(labelsText)) !== null) {
            labels.push({
                coord: [parseFloat(labelMatch[2]), parseFloat(labelMatch[1])], // [lng, lat]
                label: labelMatch[3]
            });
        }
        return labels.length > 0 ? labels : null;
    } catch (e) {
        console.warn('Failed to parse time labels:', e);
        return null;
    }
}

function hexToRgba(hex, alpha) {
    const r = parseInt(hex.slice(1, 3), 16);
    const g = parseInt(hex.slice(3, 5), 16);
    const b = parseInt(hex.slice(5, 7), 16);
    return `rgba(${r}, ${g}, ${b}, ${alpha})`;
}

function splitLineAtAntimeridian(coords) {
    // Split a line that crosses the antimeridian (±180°) into multiple segments
    // This prevents lines from going the "wrong way" around the globe
    if (!coords || coords.length < 2) return [coords];

    const segments = [];
    let currentSegment = [coords[0]];

    for (let i = 1; i < coords.length; i++) {
        const prevLng = coords[i - 1][0];
        const currLng = coords[i][0];

        // Detect antimeridian crossing (large jump in longitude)
        if (Math.abs(currLng - prevLng) > 180) {
            // Close current segment and start new one
            if (currentSegment.length > 1) {
                segments.push(currentSegment);
            }
            currentSegment = [coords[i]];
        } else {
            currentSegment.push(coords[i]);
        }
    }

    if (currentSegment.length > 1) {
        segments.push(currentSegment);
    }

    return segments;
}

function createPathPolygons(northernCoords, southernCoords) {
    // Create polygon(s) from northern and southern limits
    // Handle antimeridian crossings by creating a MultiPolygon
    if (!northernCoords || !southernCoords || northernCoords.length < 2 || southernCoords.length < 2) {
        return null;
    }

    // Split both lines at antimeridian
    const northSegments = splitLineAtAntimeridian(northernCoords);
    const southSegments = splitLineAtAntimeridian(southernCoords);

    // If no antimeridian crossing, create simple polygon
    if (northSegments.length === 1 && southSegments.length === 1) {
        const polygonCoords = [
            ...northernCoords,
            ...southernCoords.slice().reverse(),
            northernCoords[0]
        ];
        return {
            type: 'Feature',
            geometry: {
                type: 'Polygon',
                coordinates: [polygonCoords]
            }
        };
    }

    // For paths crossing antimeridian, create separate polygons for each segment pair
    // Match segments by their longitude ranges
    const polygons = [];

    for (let i = 0; i < Math.min(northSegments.length, southSegments.length); i++) {
        const northSeg = northSegments[i];
        const southSeg = southSegments[i];

        if (northSeg.length >= 2 && southSeg.length >= 2) {
            const polygonCoords = [
                ...northSeg,
                ...southSeg.slice().reverse(),
                northSeg[0]
            ];
            polygons.push([polygonCoords]);
        }
    }

    if (polygons.length === 0) return null;

    if (polygons.length === 1) {
        return {
            type: 'Feature',
            geometry: {
                type: 'Polygon',
                coordinates: polygons[0]
            }
        };
    }

    return {
        type: 'Feature',
        geometry: {
            type: 'MultiPolygon',
            coordinates: polygons
        }
    };
}

function createMultiLineString(coords) {
    // Create a LineString or MultiLineString, handling antimeridian crossings
    const segments = splitLineAtAntimeridian(coords);

    if (segments.length === 0) return null;

    if (segments.length === 1) {
        return {
            type: 'Feature',
            geometry: {
                type: 'LineString',
                coordinates: segments[0]
            }
        };
    }

    return {
        type: 'Feature',
        geometry: {
            type: 'MultiLineString',
            coordinates: segments
        }
    };
}

async function loadNasaEclipsePath(dateIso, color, prefix) {
    const eclId = formatEclipseId(dateIso);
    const nasaUrl = `https://eclipse.gsfc.nasa.gov/SEsearch/eclipse-path-data.js.php?Ecl=${encodeURIComponent(eclId)}&Spc=0.5`;
    const url = `https://dibalik.wenhop.workers.dev/?url=${encodeURIComponent(nasaUrl)}`;

    try {
        const resp = await fetch(url);
        if (!resp.ok) return false;
        const jsCode = await resp.text();

        if (jsCode.startsWith('// Error') || !jsCode.includes('centralLimitCoords')) {
            console.warn(`NASA API error for eclipse ${eclId}`);
            return false;
        }

        // Parse coordinate data from NASA's JS
        const centralCoords = parseNasaCoords(jsCode, 'centralLimitCoords');
        const northernCoords = parseNasaCoords(jsCode, 'northernLimitCoords');
        const southernCoords = parseNasaCoords(jsCode, 'southernLimitCoords');
        const timeLabels = parseNasaTimeData(jsCode);

        if (!centralCoords || centralCoords.length < 2) {
            console.warn(`No valid path data for ${eclId}`);
            return false;
        }

        // Clear any existing layers for this prefix
        clearEclipseLayers(prefix);

        // Create the totality zone polygon if we have both limits
        const pathPolygon = createPathPolygons(northernCoords, southernCoords);

        if (pathPolygon) {
            // Add polygon source
            map.addSource(`${prefix}-path-polygon`, {
                type: 'geojson',
                data: pathPolygon
            });

            // Add outer glow layer (larger, more transparent)
            map.addLayer({
                id: `${prefix}-path-fill-glow`,
                type: 'fill',
                source: `${prefix}-path-polygon`,
                paint: {
                    'fill-color': color,
                    'fill-opacity': 0.08
                }
            });

            // Add main fill layer
            map.addLayer({
                id: `${prefix}-path-fill`,
                type: 'fill',
                source: `${prefix}-path-polygon`,
                paint: {
                    'fill-color': color,
                    'fill-opacity': 0.25
                }
            });
        }

        // Add centerline with glow effect (handles antimeridian crossings)
        const centerlineGeoJson = createMultiLineString(centralCoords);
        if (centerlineGeoJson) {
            centerlineGeoJson.properties = { date: dateIso, timeLabels: timeLabels };

            map.addSource(`${prefix}-centerline`, {
                type: 'geojson',
                data: centerlineGeoJson
            });

            // Outer glow layer
            map.addLayer({
                id: `${prefix}-centerline-glow`,
                type: 'line',
                source: `${prefix}-centerline`,
                paint: {
                    'line-color': color,
                    'line-width': 12,
                    'line-opacity': 0.3,
                    'line-blur': 8
                }
            });

            // Main centerline
            map.addLayer({
                id: `${prefix}-centerline`,
                type: 'line',
                source: `${prefix}-centerline`,
                paint: {
                    'line-color': color,
                    'line-width': 3,
                    'line-opacity': 1
                }
            });
        }

        // Add limit lines as subtle dashed lines
        const northernLine = createMultiLineString(northernCoords);
        if (northernLine) {
            map.addSource(`${prefix}-northern-limit`, {
                type: 'geojson',
                data: northernLine
            });
            map.addLayer({
                id: `${prefix}-northern-limit`,
                type: 'line',
                source: `${prefix}-northern-limit`,
                paint: {
                    'line-color': color,
                    'line-width': 1.5,
                    'line-opacity': 0.6,
                    'line-dasharray': [4, 4]
                }
            });
        }

        const southernLine = createMultiLineString(southernCoords);
        if (southernLine) {
            map.addSource(`${prefix}-southern-limit`, {
                type: 'geojson',
                data: southernLine
            });
            map.addLayer({
                id: `${prefix}-southern-limit`,
                type: 'line',
                source: `${prefix}-southern-limit`,
                paint: {
                    'line-color': color,
                    'line-width': 1.5,
                    'line-opacity': 0.6,
                    'line-dasharray': [4, 4]
                }
            });
        }

        // Add hover interaction for centerline (only if centerline was added)
        if (centerlineGeoJson) {
            map.on('mouseenter', `${prefix}-centerline`, () => {
                map.getCanvas().style.cursor = 'pointer';
            });
            map.on('mouseleave', `${prefix}-centerline`, () => {
                map.getCanvas().style.cursor = '';
            });

            // Click to show time info popup
            map.on('click', `${prefix}-centerline`, (e) => {
                const coords = e.lngLat;

                // Find nearest time label
                let nearestLabel = null;
                let minDist = Infinity;

                if (timeLabels) {
                    timeLabels.forEach(item => {
                        const [lng, lat] = item.coord;
                        const dist = Math.sqrt(Math.pow(lng - coords.lng, 2) + Math.pow(lat - coords.lat, 2));
                        if (dist < minDist) {
                            minDist = dist;
                            nearestLabel = item.label;
                        }
                    });
                }

                const timeHtml = nearestLabel
                    ? `<div class="path-popup-time">${nearestLabel} UTC</div>`
                    : '';

                new mapboxgl.Popup({ closeButton: true, closeOnClick: true })
                    .setLngLat(coords)
                    .setHTML(`
                        ${timeHtml}
                        <div class="path-popup-coords">${coords.lat.toFixed(4)}°, ${coords.lng.toFixed(4)}°</div>
                        <div style="margin-top: 4px; color: ${color}; font-weight: 600;">${dateIso}</div>
                    `)
                    .addTo(map);

                e.stopPropagation();
            });
        }

        return true;
    } catch (e) {
        console.error(`Failed to load NASA path for ${dateIso}:`, e);
        return false;
    }
}

function clearPaths() {
    ['previous', 'next'].forEach(prefix => {
        clearEclipseLayers(prefix);
    });
}

async function displayResults(data) {
    let html = `<div class="coords">📍 ${data.lat.toFixed(4)}°, ${data.lon.toFixed(4)}°</div>`;

    document.getElementById('info-panel').classList.add('visible');

    html += '<div class="eclipse-card previous">';
    html += '<h3>← Previous Total Eclipse</h3>';
    if (data.previous) {
        html += `<div class="date">${data.previous.date}</div>`;
        html += `<div class="time">${data.previous.max_time_utc}</div>`;
        if (data.previous.duration_seconds) {
            html += `<div class="duration">${formatDuration(data.previous.duration_seconds)}</div>`;
        }
        html += `<a class="nasa-link" href="${getNasaUrl(data.previous.date)}" target="_blank">View on NASA</a>`;
    } else {
        html += '<p class="no-eclipse">None found in catalog</p>';
    }
    html += '</div>';

    html += '<div class="eclipse-card next">';
    html += '<h3>Next Total Eclipse →</h3>';
    if (data.next) {
        html += `<div class="date">${data.next.date}</div>`;
        html += `<div class="time">${data.next.max_time_utc}</div>`;
        if (data.next.duration_seconds) {
            html += `<div class="duration">${formatDuration(data.next.duration_seconds)}</div>`;
        }
        html += `<a class="nasa-link" href="${getNasaUrl(data.next.date)}" target="_blank">View on NASA</a>`;
    } else {
        html += '<p class="no-eclipse">None found in catalog</p>';
    }
    html += '</div>';

    document.getElementById('content').innerHTML = html;

    clearPaths();
    const legend = document.getElementById('legend');
    const loadingOverlay = document.getElementById('loading-overlay');
    let hasPath = false;

    // Show loading overlay
    loadingOverlay.style.display = 'flex';

    try {
        // Load paths in parallel
        const [prevLoaded, nextLoaded] = await Promise.all([
            data.previous ? loadNasaEclipsePath(data.previous.date, '#7eb8da', 'previous') : Promise.resolve(false),
            data.next ? loadNasaEclipsePath(data.next.date, '#f5a623', 'next') : Promise.resolve(false)
        ]);

        if (prevLoaded) hasPath = true;
        if (nextLoaded) hasPath = true;

        if (data.previous) {
            document.getElementById('legend-previous').textContent = data.previous.date;
        }
        if (data.next) {
            document.getElementById('legend-next').textContent = data.next.date;
        }
    } finally {
        // Hide loading overlay
        loadingOverlay.style.display = 'none';
    }

    legend.style.display = hasPath ? 'block' : 'none';
}
