# Sharing the LAE data with collaborators

The full archive is 48 GB: 22.3 GB of spectra and 25.6 GB of halo catalogs.
Most users only need the spectra and the halo properties of the sightlines,
and those fit in **about 22.6 GB** once the halo catalogs are trimmed to the
matched rows. The workflow is:

1. **Repackage** the tar into per-model archives plus a small matched-halo
   archive, with a checksummed `manifest.json` (`laespec pack`).
2. **Upload** the archives to a host that serves plain HTTPS downloads
   (Zenodo is recommended).
3. **Publish the manifest**: put the download URLs into `laespec/manifest.json`
   in this repo and push. After that, collaborators only need
   `pip install` and `laespec fetch`.

## 1. Repackage (on Expanse)

Run this on a compute node, or in an interactive job: it streams about 22 GB
and takes roughly 10–30 minutes on Lustre. It needs about 23 GB of free
space (about 49 GB with `--full-halos`). Use a filesystem where you have
quota headroom.

```bash
pip install --user git+https://github.com/ccain002/lae-spectra.git
laespec pack /expanse/lustre/scratch/ccain002/temp_project/for_Yongda/LAE_Catalogs.tar \
    --out /path/with/space/LAE_share \
    --full-halos            # optional: also package the 25.6 GB of full halo catalogs
```

This writes:

| archive | size |
|---|---|
| `LAE_spectra_late_end_early_start_13.tar` | 6.9 GB |
| `LAE_spectra_late_end_late_start_4.tar` | 6.9 GB |
| `LAE_spectra_late_end_very_early_start_3.tar` | 6.9 GB |
| `LAE_spectra_early_end_early_start_3.tar` | 1.6 GB |
| `LAE_matched_halos.tar` | 0.3 GB |
| `LAE_halo_catalogs_part{1,2}.tar` (with `--full-halos`) | 25.6 GB total, ≤ 20 GB each |
| `manifest.json` | sizes + SHA-256 of all of the above |

The archives are uncompressed tars. `LAEData` reads them in place, so nobody
has to extract them. Check the result with
`laespec info /path/with/space/LAE_share`.

## 2. Upload

### Option A: Zenodo (recommended)

Zenodo is free, gives you a DOI to cite in the paper, has no download
bandwidth limits, and allows up to **50 GB per record**. Spectra plus matched
halos (22.6 GB) fit in one record. Put the full halo catalogs in a second
record if you want to share them.

**Access.** Publish the record as **open** if the data can be public. That is
the simplest option, and anyone can then download it. If not, publish it with
**restricted** file access: the files are visible only to people you approve,
or to anyone holding a *share link*. Make one under the record's
*Share → Get a link* menu. For a restricted record, add the link's `token` to
the base URL in step 3, and test one download before handing it out.

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

### Option B: Globus (data stays on Expanse)

This option suits collaborators with Globus access, or the full 48 GB halo
catalogs, and avoids any upload. In the Globus web app, open the SDSC Expanse
Lustre endpoint, go to the directory with the archives, and choose
*Share → Add guest collection*. Grant read access to your collaborators'
Globus identities (or email addresses). They then transfer the files to their
own endpoint. `laespec fetch` doesn't apply here: they point `LAEData` at
wherever the files land.

Lustre **scratch is purged**. Before using this option, move the data to
project storage or a location that won't be deleted.

### Option C: collaborators who already have Expanse accounts

Make the tar world-readable (it already is: `-rwxrwxrwx`) and every parent
directory traversable (`chmod o+x` on each). Collaborators then run
`LAEData("/expanse/.../LAE_Catalogs.tar")` directly, with no copying.

## 3. Publish the manifest

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
laespec fetch --dest /tmp/lae_test --model early_end_early_start_3   # smallest archive
laespec info /tmp/lae_test
```

## 4. Give collaborators access

- **Code**: the repo is private. Either add each person
  (`gh api -X PUT repos/ccain002/lae-spectra/collaborators/<github-username>`,
  or *Settings → Collaborators* on GitHub), or make it public with
  `gh repo edit ccain002/lae-spectra --visibility public --accept-visibility-change-consequences`.
  For a private repo, collaborators install with the same `pip install git+https://…` command
  after authenticating to GitHub.
- **Data**: send them the Zenodo record or share link, or add them to the
  Globus collection.
- **Instructions for them**: the *Install*, *Get the data* and *Usage*
  sections of the README.

Message you can send:

> Code + docs: https://github.com/ccain002/lae-spectra
> ```
> pip install git+https://github.com/ccain002/lae-spectra.git
> laespec fetch --dest ~/lae_data            # ~22.6 GB; add --model NAME for a single model
> python -c "import laespec as ls; print(ls.LAEData('~/lae_data').summary())"
> ```
> See the README for the data format and examples.

## Updating the data later

Rerun `laespec pack`. On Zenodo, create a **new version** of the record: it
gets a new record ID while the concept DOI stays the same. Then rerun
`set-urls` with the new ID and push the updated manifest. Old installs keep
pointing at the old version until people `pip install -U`.
