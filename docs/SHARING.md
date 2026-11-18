# Sharing the LAE data with collaborators

The data is shared in two parts:

- **The bundled example in the repo.** One real spectrum file
  (`late_end_early_start_13`, z = 6.9746, `skip=539`, 40 MB) is in
  `laespec/example_data/`, along with the wavelength grid and the matched
  halo rows. It comes with every clone and every `pip install`, so
  `laespec.example()`, the quickstart and the tests work right away.
- **The full 45.5 GB dataset on Globus.** That's 20.4 GB of spectra plus
  25.1 GB of halo catalogs.

## Globus

The guest collection is live:
<https://app.globus.org/file-manager?origin_id=e06caa06-533b-47fb-b874-7c120d8cf896&origin_path=%2F>
(collection ID `e06caa06-533b-47fb-b874-7c120d8cf896`, rooted at
`/expanse/lustre/projects/ucr138/ccain002/LAE_Catalogs/`). To change who can
see it, use the collection's **Permissions** tab in the Globus web app.

How it was set up:

1. **The data lives on project storage**, not scratch, because Lustre scratch is purged:
   `/expanse/lustre/projects/ucr138/ccain002/LAE_Catalogs/LAE_Catalogs.tar`.
   That directory holds only this tar, so a guest collection rooted there
   shares nothing else. `LAEData` reads the tar as-is. Repackaging
   (`laespec pack`, below) is optional: it only helps people who want just
   one model (about 7 GB) instead of the whole 45.5 GB.
2. In the Globus web app, open the SDSC Expanse endpoint and go to
   `/expanse/lustre/projects/ucr138/ccain002/LAE_Catalogs/`. Choose *Share → Add guest collection*.
3. Grant **read** access to your collaborators' Globus identities (or email
   addresses).
4. Collaborators transfer the tar to their own endpoint (a Globus Connect
   Personal endpoint works on a laptop), then run
   `LAEData("/where/it/landed/LAE_Catalogs.tar")`.

## Collaborators with Expanse accounts

The tar is world-readable and every parent directory is traversable, so any
Expanse user can read it in place with no copying:
`LAEData("/expanse/lustre/projects/ucr138/ccain002/LAE_Catalogs/LAE_Catalogs.tar")`.

## Code access

The repo is public, so anyone can clone it or `pip install git+https://…` it
without a GitHub account. Only the full dataset needs an explicit grant on
Globus.

Message you can send:

> Code + docs: https://github.com/ccain002/lae-spectra
> ```
> pip install git+https://github.com/ccain002/lae-spectra.git   # includes a 40 MB example file
> python -c "import laespec as ls; print(ls.example().spectra('late_end_early_start_13', 7.0).flux.shape)"
> ```
> The full dataset (45.5 GB) is on Globus:
> https://app.globus.org/file-manager?origin_id=e06caa06-533b-47fb-b874-7c120d8cf896&origin_path=%2F
> Transfer LAE_Catalogs.tar and use `ls.LAEData("/path/to/LAE_Catalogs.tar")`.
> See the README for the data format and examples.

## Alternative: Zenodo + `laespec fetch`

Use this route if you later want a public, citable copy with a DOI, or plain
HTTPS downloads that need no Globus. `laespec pack` splits the tar into
per-model archives with a checksummed manifest. `laespec fetch` then
downloads them, with resume support. This route isn't set up yet:
`laespec/manifest.json` is still empty, so `laespec fetch` has nothing to
download.

### A1. Repackage (on Expanse)

Run this on a compute node, or in an interactive job: it streams about 21 GB
and takes roughly 10–30 minutes on Lustre. It needs about 21 GB of free
space (about 46 GB with `--full-halos`). Use a filesystem where you have
quota headroom.

```bash
pip install --user git+https://github.com/ccain002/lae-spectra.git
laespec pack /expanse/lustre/projects/ucr138/ccain002/LAE_Catalogs/LAE_Catalogs.tar \
    --out /path/with/space/LAE_share \
    --full-halos            # optional: also package the 25.1 GB of full halo catalogs
```

This writes:

| archive | size |
|---|---|
| `LAE_spectra_late_end_early_start_13.tar` | 6.8 GB |
| `LAE_spectra_late_end_late_start_4.tar` | 6.8 GB |
| `LAE_spectra_late_end_very_early_start_3.tar` | 6.8 GB |
| `LAE_matched_halos.tar` | 0.3 GB |
| `LAE_halo_catalogs_part{1,2}.tar` (with `--full-halos`) | 25.1 GB total, ≤ 20 GB each |
| `manifest.json` | sizes + SHA-256 of all of the above |

The archives are uncompressed tars. `LAEData` reads them in place, so nobody
has to extract them. Check the result with
`laespec info /path/with/space/LAE_share`.

### A2. Upload to Zenodo

Zenodo is free, gives you a DOI to cite in the paper, has no download
bandwidth limits, and allows up to **50 GB per record**. Spectra plus matched
halos (20.7 GB) fit in one record. Put the full halo catalogs in a second
record if you want to share them.

**Access.** Publish the record as **open** if the data can be public. That is
the simplest option, and anyone can then download it. If not, publish it with
**restricted** file access: the files are visible only to people you approve,
or to anyone holding a *share link*. Make one under the record's
*Share → Get a link* menu. For a restricted record, add the link's `token` to
the base URL in step A3, and test one download before handing it out.

Upload straight from Expanse with the REST API, so the files don't have to
pass through your laptop:

```bash
# 1. Create a personal token at https://zenodo.org/account/settings/applications
#    with scopes deposit:write and deposit:actions
export ZTOKEN=...            # do not commit this anywhere

# 2. Create an empty deposition; note "id" and "links.bucket" in the reply
curl -s -H "Authorization: Bearer $ZTOKEN" -H "Content-Type: application/json" \
     -X POST https://zenodo.org/api/deposit/depositions -d '{}' | python -m json.tool

export BUCKET=https://zenodo.org/api/files/<bucket-id>   # from links.bucket
export RECORD=<id>

# 3. Upload each archive (each takes a few minutes; run them in a screen/tmux session)
cd /path/with/space/LAE_share
for f in LAE_spectra_*.tar LAE_matched_halos.tar; do
  curl --fail --progress-bar -o /dev/null -H "Authorization: Bearer $ZTOKEN" \
       --upload-file "$f" "$BUCKET/$f"
done
```

Then open `https://zenodo.org/deposit/$RECORD` in a browser. Fill in the
title, authors, description (copy the *Data format* section of the README),
license and access level, and press **Publish**.

To do a dry run first, use `https://sandbox.zenodo.org` everywhere above. It
needs a separate account and token.

### A3. Publish the manifest

Fill in the URLs and commit the manifest into the package, so every install
carries it:

```bash
# Zenodo, open access:
laespec set-urls /path/with/space/LAE_share/manifest.json \
    --base-url "https://zenodo.org/records/$RECORD/files/{name}?download=1"
# Zenodo, restricted access with a share link:
#   --base-url "https://zenodo.org/records/$RECORD/files/{name}?download=1&token=<share-token>"

cp /path/with/space/LAE_share/manifest.json laespec/manifest.json
git add laespec/manifest.json && git commit -m "Add data manifest" && git push
```

If you uploaded the halo catalogs to a different record, run `set-urls` on a
copy for that record and merge the `files` entries by hand. Each entry carries
its own `url`.

Test it as a collaborator would:

```bash
pip install --force-reinstall git+https://github.com/ccain002/lae-spectra.git
laespec fetch --list --dest /tmp/lae_test
laespec fetch --dest /tmp/lae_test --model late_end_late_start_4   # one model, 6.8 GB
laespec info /tmp/lae_test
```

### Updating the data later

Rerun `laespec pack`. On Zenodo, create a **new version** of the record: it
gets a new record ID while the concept DOI stays the same. Then rerun
`set-urls` with the new ID and push the updated manifest. Old installs keep
pointing at the old version until people `pip install -U`.
