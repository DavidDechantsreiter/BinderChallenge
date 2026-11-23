# BinderChallenge

Commands are for Apple Silicon Macs; Windows and Linux need different kubectl and runai downloads.
## Step 1: Access.
- Get your GASPAR username from IS-Academia. 
- Check onto the jump host with id `ssh <gaspar>@jumphost.rcp.epfl.ch`. 
- Once connected type `id`. Look for `rcp-caas-make-dlai-g11_AppGrpU` and `rcp-runai-make-dlai_AppGrpU`.
- Additionally, note your `uid` and `gid` it's important for point 5 of setting up

## Step 2: VPN.
- In order to work also remotely, you need to connect to the EPFL network via VPN. EPFL uses Cisco Secure Client (https://www.epfl.ch/campus/services/en/it-services/network-services/remote-intranet-access/vpn-clients-available/).
- Once set up connect to `vpn.epfl.ch`, login using <gaspar>@epfl.ch + GASPAR password + Authenticator code (no notification, just the code on the screen)

You now have access to the jumphost. The jumphost is not meant for heavy coding or modelling work. Instead it serves as a quick way to check and move files around, clone repositories and download data. Next you will see how to set up the actual access to the cluster so that you can use the gpu and cpus required for the project.

## Step 3: What is RCP
EPFL's Research Computing Platform (RCP) is a shared cluster: a large pool of servers with GPUs of different generations (V100, A100, H100, ...). It's managed by **Kubernetes** and **Run:ai**, which decide who gets which GPU and bill the time to our project's credits.

# Step 4: Set up runai and kubernetes
1. Install the tools (Mac terminal):
```bash
    curl -LO "https://dl.k8s.io/release/v1.30.11/bin/darwin/arm64/kubectl"
    chmod +x kubectl && sudo mv kubectl /usr/local/bin/kubectl
    mkdir -p ~/.kube
    curl -o ~/.kube/config https://raw.githubusercontent.com/epfml/getting-started/main/kubeconfig.yaml
    curl -L -o runai https://rcp-caas-prod.rcp.epfl.ch/cli/darwin
    chmod +x runai && sudo mv runai /usr/local/bin/runai
```
    sudo wants the Mac password. If macOS blocks runai: xattr -d com.apple.quarantine /usr/local/bin/runai.

2. Log in and pick the project:
```bash
    runai login                          # browser, EPFL account
    runai list projects
    runai config project <dlai-project>  # DLAI as default (fine if it's your only project, Jacob don't run this)
```

To run anything on the cluster, you need an **image**: a ready-made Linux system (operating system, CUDA, Python, tools and your user account) packed into one file. Every time you start a job, the server launches a fresh copy of your image, called a **container**, and that's the environment you work in.

You make an image in two steps:
1. **Write a Dockerfile**: a text file of setup instructions ("start from Ubuntu with CUDA, install git and Python, create my user"), I've already created this for you.
2. **Build it** on your laptop: Docker follows the instructions and saves the result as an image, which you upload to EPFL's registry (basically a database of images, so the cluster can download it from there).

You only build once. It's slow (20+ minutes) because it installs everything, but afterwards every job starts from the finished image in a minute or two. The models will require their own Python packages, but we'll handle that once were on the cluster and make a virtual enviornment that lives in cluster storage instead.

3. Build your image
- Install Docker Desktop (this will help you create and manage your image): ```Docker Desktop: brew install --cask docker``` open it once. Remember the `gid` and `uid` you got from running `id` on the jumphost. Navigate to the repo's Dockerfile folder:
    ```bash
    docker login registry.rcp.epfl.ch
    docker buildx build --platform linux/amd64 \
        --build-arg NB_USER=<gaspar> --build-arg NB_UID=<uid> \
        --build-arg NB_GROUP=<dlai-group> --build-arg NB_GID=<dlai-gid> \
        -t registry.rcp.epfl.ch/<dlai-registry>/<gaspar>:v1 --push .
    ```
    --platform linux/amd64 is required because Macs are ARM and the cluster is x86. Expect 20+ minutes the first time.
7.  Test pod: submit (see Daily workflow), then runai exec ... -- bash, id (their name and UID), nvidia-smi, touch a file in their folder.
8.  GitHub SSH key on the jump host: ssh-keygen -t ed25519, add ~/.ssh/id_ed25519.pub to GitHub, clone with the git@github.com: URL.
9.  venv + init script on DLAI storage, as in step 5–6 of your steps.
For the Dockerfile: in the team repo, replace your hard-coded user lines with build arguments at the top, so step 5 works for everyone:
ARG NB_USER NB_UID NB_GROUP NB_GID
No passwords or tokens go in the Dockerfile or the image.

You don't log in to a GPU machine directly. Instead:

1. **You submit a job** from your laptop (`runai submit`), saying how many GPUs/CPUs you need and which **image** to use.
- The image is a frozen, ready-made Linux system: OS, CUDA, Python and tools are already inside it, so nothing gets installed at startup.
2. **The cluster finds a free server** and starts a **container** from your image. 
3. **You connect to the container** (`runai exec`, or VS Code) and work there as if it were a remote Linux machine.

**The container is temporary.** When the job ends, everything inside it is wiped. Only files on the **mounted storage** (our team folder, e.g. `/dlai`) survive. So keep code, data, venvs and results there, never just "in the container".

**Delete your job when you're done.** A running job keeps its GPU reserved and keeps billing credits, even when idle.