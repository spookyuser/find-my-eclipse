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

    document.getElementById('summary-tab').innerHTML = '<p class="loading">Loading...</p>';

    // Reset to Summary tab and fetch all eclipses in parallel
    switchToTab('summary');
    fetchAllEclipses(lat, lng);

    fetch(`/api/eclipses?lat=${lat.toFixed(6)}&lon=${lng.toFixed(6)}`)
        .then(r => r.json())
        .then(data => {
            displayResults(data);
        })
        .catch(err => {
            document.getElementById('summary-tab').innerHTML = `<p class="no-eclipse">Error: ${err.message}</p>`;
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
    // Convert ISO date to NASA format
    // CE dates: 2024-04-08 -> +20240408
    // BCE dates: -0549-06-30 -> -05490630
    if (dateIso.startsWith('-')) {
        // BCE date - preserve the minus sign, remove only the date separators
        return '-' + dateIso.slice(1).replace(/-/g, '');
    }
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
    const url = `https://eclipse.wenhop.workers.dev/?url=${encodeURIComponent(nasaUrl)}`;

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
                        <div style="margin-top: 4px; color: ${color}; font-weight: 600;">${formatDisplayDate(dateIso)}</div>
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
        html += `<div class="date">${formatDisplayDate(data.previous.date)}</div>`;
        html += `<a class="nasa-link" href="${getNasaUrl(data.previous.date)}" target="_blank">View on NASA</a>`;
    } else {
        html += '<p class="no-eclipse">None found in catalog</p>';
    }
    html += '</div>';

    html += '<div class="eclipse-card next">';
    html += '<h3>Next Total Eclipse →</h3>';
    if (data.next) {
        html += `<div class="date">${formatDisplayDate(data.next.date)}</div>`;
        html += `<a class="nasa-link" href="${getNasaUrl(data.next.date)}" target="_blank">View on NASA</a>`;
    } else {
        html += '<p class="no-eclipse">None found in catalog</p>';
    }
    html += '</div>';

    document.getElementById('summary-tab').innerHTML = html;

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
            document.getElementById('legend-previous').textContent = formatDisplayDate(data.previous.date);
        }
        if (data.next) {
            document.getElementById('legend-next').textContent = formatDisplayDate(data.next.date);
        }
    } finally {
        // Hide loading overlay
        loadingOverlay.style.display = 'none';
    }

    legend.style.display = hasPath ? 'block' : 'none';
}

// Tab switching
function switchToTab(tabName) {
    document.querySelectorAll('.tab-btn').forEach(btn => {
        btn.classList.toggle('active', btn.dataset.tab === tabName);
    });
    document.querySelectorAll('.tab-panel').forEach(panel => {
        panel.classList.toggle('active', panel.id === `${tabName}-tab`);
    });
}

document.querySelectorAll('.tab-btn').forEach(btn => {
    btn.addEventListener('click', () => {
        switchToTab(btn.dataset.tab);
    });
});

// All eclipses functionality
let allEclipsesData = null;
let currentSelectedEclipse = null;

async function fetchAllEclipses(lat, lng) {
    const listEl = document.getElementById('eclipse-list');
    const loadingEl = document.getElementById('eclipse-list-loading');
    const chartLoadingEl = document.getElementById('chart-loading');
    const chartEl = document.getElementById('gap-chart');

    listEl.innerHTML = '';
    chartEl.innerHTML = '';
    loadingEl.style.display = 'block';
    chartLoadingEl.style.display = 'block';
    allEclipsesData = null;

    try {
        const resp = await fetch(`/api/all-eclipses?lat=${lat.toFixed(6)}&lon=${lng.toFixed(6)}`);
        const data = await resp.json();
        allEclipsesData = data;
        renderEclipseList(data);
        renderGapChart(data);
    } catch (err) {
        listEl.innerHTML = `<p class="no-eclipse">Error: ${err.message}</p>`;
        chartEl.innerHTML = `<p class="no-eclipse">Error: ${err.message}</p>`;
        chartLoadingEl.style.display = 'none';
    } finally {
        loadingEl.style.display = 'none';
    }
}

function formatDisplayDate(dateIso) {
    const months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];

    if (dateIso.startsWith('-')) {
        // BCE date: -0549-06-30 -> Jun 30, 549 BCE
        const parts = dateIso.slice(1).split('-');
        const year = parseInt(parts[0], 10);
        const month = months[parseInt(parts[1], 10) - 1];
        const day = parseInt(parts[2], 10);
        return `${month} ${day}, ${year} BCE`;
    }

    // CE date: 2024-04-08 -> Apr 8, 2024
    const parts = dateIso.split('-');
    const year = parseInt(parts[0], 10);
    const month = months[parseInt(parts[1], 10) - 1];
    const day = parseInt(parts[2], 10);
    return `${month} ${day}, ${year}`;
}

function getYearsAgo(dateIso) {
    const currentYear = new Date().getFullYear();
    let eclipseYear;

    if (dateIso.startsWith('-')) {
        // BCE year (negative)
        eclipseYear = -parseInt(dateIso.slice(1).split('-')[0], 10);
    } else {
        eclipseYear = parseInt(dateIso.split('-')[0], 10);
    }

    const diff = currentYear - eclipseYear;
    if (diff > 0) {
        return `${diff.toLocaleString()} years ago`;
    } else if (diff < 0) {
        return `in ${Math.abs(diff).toLocaleString()} years`;
    }
    return 'this year';
}

function getYearFromDate(dateIso) {
    if (dateIso.startsWith('-')) {
        return -parseInt(dateIso.slice(1).split('-')[0], 10);
    }
    return parseInt(dateIso.split('-')[0], 10);
}

function renderGapChart(data) {
    const chartEl = document.getElementById('gap-chart');
    const loadingEl = document.getElementById('chart-loading');

    loadingEl.style.display = 'none';

    if (!data || !data.eclipses || data.eclipses.length < 2) {
        chartEl.innerHTML = '<p class="no-eclipse">Not enough eclipses to show gaps</p>';
        return;
    }

    const eclipses = data.eclipses;
    const currentYear = new Date().getFullYear();

    // Calculate gaps between consecutive eclipses
    const gaps = [];
    for (let i = 0; i < eclipses.length - 1; i++) {
        const year1 = getYearFromDate(eclipses[i].date);
        const year2 = getYearFromDate(eclipses[i + 1].date);
        gaps.push({
            startYear: year1,
            endYear: year2,
            startDate: eclipses[i].date,
            endDate: eclipses[i + 1].date,
            gapYears: year2 - year1
        });
    }

    // Get min and max years for the timeline
    const minYear = getYearFromDate(eclipses[0].date);
    const maxYear = getYearFromDate(eclipses[eclipses.length - 1].date);
    const yearRange = maxYear - minYear;

    // SVG dimensions
    const width = chartEl.clientWidth || 320;
    const height = 200;
    const padding = { top: 40, right: 20, bottom: 30, left: 20 };
    const chartWidth = width - padding.left - padding.right;
    const chartHeight = height - padding.top - padding.bottom;

    // Scale function to convert year to x position
    const yearToX = (year) => {
        return padding.left + ((year - minYear) / yearRange) * chartWidth;
    };

    // Build SVG
    let svg = `<svg class="gap-chart-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet">`;

    // Draw axis line
    svg += `<line class="chart-axis" x1="${padding.left}" y1="${padding.top + chartHeight/2}" x2="${width - padding.right}" y2="${padding.top + chartHeight/2}" />`;

    // Draw gaps as colored bars
    const barHeight = 24;
    const barY = padding.top + chartHeight/2 - barHeight/2;

    // Collect gradient definitions for gaps spanning the current year
    const gradientDefs = [];

    gaps.forEach((gap, i) => {
        const x1 = yearToX(gap.startYear);
        const x2 = yearToX(gap.endYear);
        const isPast = gap.endYear <= currentYear;
        const isFuture = gap.startYear >= currentYear;
        const isCurrent = gap.startYear < currentYear && gap.endYear > currentYear;

        let color;
        if (isPast) {
            color = '#7eb8da';
        } else if (isFuture) {
            color = '#f5a623';
        } else {
            // Gap spans current year - calculate where "Now" falls within the gap
            const nowPercent = ((currentYear - gap.startYear) / (gap.endYear - gap.startYear)) * 100;
            const gradientId = `gradientMixed${i}`;
            gradientDefs.push(`
                <linearGradient id="${gradientId}" x1="0%" y1="0%" x2="100%" y2="0%">
                    <stop offset="0%" style="stop-color:#7eb8da"/>
                    <stop offset="${nowPercent}%" style="stop-color:#7eb8da"/>
                    <stop offset="${nowPercent}%" style="stop-color:#f5a623"/>
                    <stop offset="100%" style="stop-color:#f5a623"/>
                </linearGradient>
            `);
            color = `url(#${gradientId})`;
        }

        const gapWidth = Math.max(x2 - x1, 2);

        svg += `<rect class="gap-bar"
            x="${x1}" y="${barY}"
            width="${gapWidth}" height="${barHeight}"
            fill="${color}"
            rx="3"
            data-gap-index="${i}"
            data-years="${gap.gapYears}"
            data-start="${gap.startDate}"
            data-end="${gap.endDate}" />`;

        // Add year label inside the bar if wide enough
        const labelText = `${gap.gapYears}y`;
        const minWidthForLabel = labelText.length * 7 + 10; // Approximate text width
        if (gapWidth >= minWidthForLabel) {
            const labelX = x1 + gapWidth / 2;
            const labelY = barY + barHeight / 2 + 4;
            svg += `<text class="gap-bar-label" x="${labelX}" y="${labelY}" text-anchor="middle" fill="#fff" font-size="10" font-weight="500" pointer-events="none">${gap.gapYears}y</text>`;
        }
    });

    // Add gradient definitions for mixed past/future gaps
    if (gradientDefs.length > 0) {
        svg += `<defs>${gradientDefs.join('')}</defs>`;
    }

    // Draw eclipse markers on top
    eclipses.forEach((ecl, i) => {
        const year = getYearFromDate(ecl.date);
        const x = yearToX(year);
        const isPast = year < currentYear;
        const color = isPast ? '#7eb8da' : '#f5a623';

        svg += `<circle class="eclipse-marker"
            cx="${x}" cy="${padding.top + chartHeight/2}" r="4"
            fill="${color}" stroke="#fff" stroke-width="1" />`;
    });

    // Draw "now" line if within range
    if (currentYear >= minYear && currentYear <= maxYear) {
        const nowX = yearToX(currentYear);
        svg += `<line class="now-line" x1="${nowX}" y1="${padding.top}" x2="${nowX}" y2="${padding.top + chartHeight}" />`;
        svg += `<text x="${nowX}" y="${padding.top - 8}" fill="#fff" font-size="10" text-anchor="middle">Now</text>`;
    }

    // Add year labels at start and end
    const formatYearLabel = (year) => {
        if (year < 0) {
            return `${Math.abs(year)} BCE`;
        }
        return year.toString();
    };

    svg += `<text class="gap-label" x="${padding.left}" y="${height - 8}" text-anchor="start">${formatYearLabel(minYear)}</text>`;
    svg += `<text class="gap-label" x="${width - padding.right}" y="${height - 8}" text-anchor="end">${formatYearLabel(maxYear)}</text>`;

    // Add some intermediate labels if space allows
    if (chartWidth > 200) {
        const midYear = Math.round((minYear + maxYear) / 2);
        svg += `<text class="gap-label" x="${yearToX(midYear)}" y="${height - 8}" text-anchor="middle">${formatYearLabel(midYear)}</text>`;
    }

    svg += '</svg>';

    // Add header with data range note
    chartEl.innerHTML = `<div class="chart-header">Total Solar Eclipses<span class="chart-subheader">${eclipses.length} eclipses • 2000 BCE – 3000 CE</span></div>${svg}`;

    // Add tooltip container
    const tooltip = document.createElement('div');
    tooltip.className = 'gap-tooltip';
    tooltip.style.display = 'none';
    chartEl.style.position = 'relative';
    chartEl.appendChild(tooltip);

    // Add hover interactions for gaps
    chartEl.querySelectorAll('.gap-bar').forEach(bar => {
        bar.addEventListener('mouseenter', (e) => {
            const years = bar.dataset.years;
            const startDate = bar.dataset.start;
            const endDate = bar.dataset.end;

            tooltip.innerHTML = `
                <div class="gap-tooltip-years">${years} year${years !== '1' ? 's' : ''} between eclipses</div>
                <div class="gap-tooltip-dates">${formatDisplayDate(startDate)} → ${formatDisplayDate(endDate)}</div>
            `;
            tooltip.style.display = 'block';
        });

        bar.addEventListener('mousemove', (e) => {
            const rect = chartEl.getBoundingClientRect();
            let left = e.clientX - rect.left + 10;
            let top = e.clientY - rect.top - 40;

            // Keep tooltip in bounds
            const tooltipRect = tooltip.getBoundingClientRect();
            if (left + tooltipRect.width > rect.width) {
                left = e.clientX - rect.left - tooltipRect.width - 10;
            }
            if (top < 0) {
                top = e.clientY - rect.top + 20;
            }

            tooltip.style.left = left + 'px';
            tooltip.style.top = top + 'px';
        });

        bar.addEventListener('mouseleave', () => {
            tooltip.style.display = 'none';
        });
    });
}

function renderEclipseList(data) {
    const listEl = document.getElementById('eclipse-list');

    if (!data.eclipses || data.eclipses.length === 0) {
        listEl.innerHTML = '<p class="no-eclipse">No total eclipses found at this location</p>';
        return;
    }

    const now = new Date();
    let html = `<div class="eclipse-list-header">Total Solar Eclipses<span class="chart-subheader">${data.total_count} eclipses • 2000 BCE – 3000 CE</span></div>`;
    html += '<div class="eclipse-list-items">';

    data.eclipses.forEach(ecl => {
        const isPast = ecl.date.startsWith('-') || new Date(ecl.date) < now;
        const durationText = ecl.duration_seconds
            ? `${Math.floor(ecl.duration_seconds / 60)}m ${Math.round(ecl.duration_seconds % 60)}s`
            : '';

        html += `
            <div class="eclipse-list-item ${isPast ? 'past' : 'future'}"
                 data-cat-no="${ecl.cat_no}"
                 data-date="${ecl.date}"
                 title="${getYearsAgo(ecl.date)}">
                <div class="eclipse-list-item-date">${formatDisplayDate(ecl.date)}</div>
                ${durationText ? `<div class="eclipse-list-item-duration">${durationText}</div>` : ''}
                <div class="eclipse-list-item-badge">${isPast ? 'Past' : 'Future'}</div>
            </div>
        `;
    });

    html += '</div>';
    listEl.innerHTML = html;

    // Add click handlers
    listEl.querySelectorAll('.eclipse-list-item').forEach(item => {
        item.addEventListener('click', () => {
            const catNo = parseInt(item.dataset.catNo, 10);
            const date = item.dataset.date;
            showEclipseFromList(catNo, date, item);
        });
    });
}

async function showEclipseFromList(catNo, date, itemEl) {
    // Update selection visual
    document.querySelectorAll('.eclipse-list-item').forEach(el => {
        el.classList.remove('selected');
    });
    itemEl.classList.add('selected');

    // Determine if eclipse is past or future
    const isPast = date.startsWith('-') || new Date(date) < new Date();
    const color = isPast ? '#7eb8da' : '#f5a623';
    const prefix = isPast ? 'previous' : 'next';

    // Clear existing paths and show loading
    clearPaths();
    const loadingOverlay = document.getElementById('loading-overlay');
    loadingOverlay.style.display = 'flex';

    try {
        // Load the selected eclipse path
        await loadNasaEclipsePath(date, color, prefix);

        // Update legend
        const legend = document.getElementById('legend');
        const previousLegend = document.getElementById('legend-previous').parentElement;
        const nextLegend = document.getElementById('legend-next').parentElement;

        if (isPast) {
            document.getElementById('legend-previous').textContent = formatDisplayDate(date);
            previousLegend.style.display = '';
            nextLegend.style.display = 'none';
        } else {
            document.getElementById('legend-next').textContent = formatDisplayDate(date);
            nextLegend.style.display = '';
            previousLegend.style.display = 'none';
        }
        legend.style.display = 'block';
    } finally {
        loadingOverlay.style.display = 'none';
    }
}
