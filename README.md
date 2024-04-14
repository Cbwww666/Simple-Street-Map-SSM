# Simple-Street-Map-SSM-
simple terminal-based street map application named Simple Street Map (SSM)

## Background
OpenStreetMapLinks to an external site. (OSM) is an open source mapping platform that allows people to contribute, edit, and use geospatial data freely. The platform employs user-generated data to create detailed and up-to-date maps, encompassing a diverse range of geographical features such as roads, buildings, rivers, and points of interest. OSM has gained prominence for its flexibility, inclusivity, and the ability to serve various applications, from navigation tools and urban planning to disaster response and humanitarian aid efforts.

## Simple Street Map
There are only two main concepts in simple street map: nodes and ways.

### Node
A node in OSM is a single point on the map, denoted by latitude and longitude. A node can optionally contain metadata, but we have removed all of them for the simplicity of the assignment, e.g., turn restrictions like "no left turn".

### Way
A way in OSM contains an ordered sequence of nodes that makes up a road segment. It also contains various metadata about the road segment. For this assignment, we only care about three pieces of information: 

1. Name: the name of the road segment. Note that some roads are unnamed. In this case, they are given a unique OSM way ID.
2. One-way: whether the road segment is one-way, signifying that you cannot drive in the reverse direction.
3. Speed Limit: the maximum legal speed of the road segment.

### Example
![image](https://github.com/Cbwww666/Simple-Street-Map-SSM/assets/67548133/e335c0f5-a0e6-452d-8463-590e46a2c87c)

