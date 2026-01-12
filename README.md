# Find My Eclipse

To be honest i have no idea how this works and i especially don't understand how the stuff in [eclipse/compute.py](eclipse/compute.py) (which gpt-5.2 pro wrote) works but it does seem to work.

If you have any ideas on how to improve it, please let me know!

# Run

```bash
uv sync
uvicorn eclipse.web:app --reload
```

# View

[http://localhost:8000](http://localhost:8000)


# Command Line  

```bash
uv run eclipse --lat 51.5074 --lon -0.1278
```

# Acknowledgments

- [NASA's GSFC, Fred Espenak](https://eclipse.gsfc.nasa.gov/) for the eclipse data
- gpt-5.2-pro for writing 80% of the code
- claude-4.5-opus for writing 20% of the code
- The sun, moon and earth 

