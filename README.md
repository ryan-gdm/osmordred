## **Osmordred**: Unified RDKit Descriptors in C++

### Build and Start the Container

```bash
docker-compose up -d --build
```

### Build the Extension

Enter the container:

```bash
docker exec -it osmordred bash
```

Then inside the container:

```bash
bash build.sh
```

### Run Tests

```bash
cd test
conda activate osmordred
python test.py
```
