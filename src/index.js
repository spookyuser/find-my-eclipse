// src/index.js

const CORS_MAX_AGE = "86400";

function teapot() {
  return new Response(null, { status: 418, statusText: "I'm a teapot" });
}

function isPreflight(request) {
  const h = request.headers;
  return (
    h.get("Origin") !== null &&
    h.get("Access-Control-Request-Method") !== null
  );
}

function handleOptions(request) {
  if (!isPreflight(request)) {
    return new Response(null, {
      headers: {
        Allow: "OPTIONS, GET, HEAD, POST, PUT, PATCH, DELETE",
      },
    });
  }

  const origin = request.headers.get("Origin");
  const reqHeaders = request.headers.get("Access-Control-Request-Headers") || "";
  const reqMethod = request.headers.get("Access-Control-Request-Method") || "";

  return new Response(null, {
    status: 204,
    headers: {
      "Access-Control-Allow-Origin": origin,
      "Access-Control-Allow-Credentials": "true",
      "Access-Control-Allow-Methods":
        reqMethod || "GET,HEAD,POST,PUT,PATCH,DELETE,OPTIONS",
      "Access-Control-Allow-Headers": reqHeaders,
      "Access-Control-Max-Age": CORS_MAX_AGE,
      Vary: "Origin",
    },
  });
}

function addCorsHeaders(response, origin) {
  response = new Response(response.body, response);
  if (origin !== null) {
    response.headers.set("Access-Control-Allow-Origin", origin);
    response.headers.set("Access-Control-Allow-Credentials", "true");
    response.headers.set(
      "Access-Control-Allow-Methods",
      "GET,HEAD,POST,PUT,PATCH,DELETE,OPTIONS"
    );
    response.headers.append("Vary", "Origin");
  }
  return response;
}

async function fetchFromUpstream(request, target) {
  const upstreamUrl = new URL(target);
  const upstreamReq = new Request(upstreamUrl.toString(), request);
  upstreamReq.headers.set("Origin", upstreamUrl.origin);
  return fetch(upstreamReq);
}

async function handleRequest(request, ctx) {
  const origin = request.headers.get("Origin");
  const url = new URL(request.url);
  const target = url.searchParams.get("url");
  if (target == null) return teapot();

  // Only cache GET/HEAD requests
  if (request.method !== "GET" && request.method !== "HEAD") {
    const response = await fetchFromUpstream(request, target);
    return addCorsHeaders(response, origin);
  }

  const cache = caches.default;

  // Use the full worker URL as cache key (includes ?url=... param)
  const cacheKey = new Request(url.toString(), { method: "GET" });

  // Check cache first
  let response = await cache.match(cacheKey);

  if (!response) {
    // Cache miss - fetch from upstream
    response = await fetchFromUpstream(request, target);

    // Only cache successful responses
    if (response.ok) {
      // Clone response (body can only be read once)
      const responseToCache = response.clone();

      // Set long cache TTL (1 year)
      const cachedResponse = new Response(responseToCache.body, responseToCache);
      cachedResponse.headers.set("Cache-Control", "public, max-age=31536000, immutable");

      // Store in cache (don't await - fire and forget)
      ctx.waitUntil(cache.put(cacheKey, cachedResponse));
    }
  }

  return addCorsHeaders(response, origin);
}

export default {
  async fetch(request, env, ctx) {
    if (request.method === "OPTIONS") return handleOptions(request);
    return handleRequest(request, ctx);
  },
};
