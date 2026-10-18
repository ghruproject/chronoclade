# One-time Pasteur authorisation

Native Klebsiella lineage typing needs a reference database from BIGSdb-Pasteur.
A client key identifies your application; it does not replace your account's
approval to access the database. `authorize-pasteur` handles this approval through
your own browser and saves a private access token for the installed typing tools.
It does not sign in, submit a consent form or open a browser automatically.

## Install, approve, then build

Install the tools first, using the destination you will use throughout:

```bash
pixi run chronoclade setup-typing --species klebsiella \
  --destination ~/.local/share/chronoclade/typing
```

Keep the issued Pasteur client-key letter or supported JSON in a private local
file. Supply its **path**, rather than pasting keys into the command line:

```bash
pixi run chronoclade authorize-pasteur \
  --destination ~/.local/share/chronoclade/typing \
  --pasteur-secrets /path/to/pasteur-client-key.txt
```

Open the printed `https://bigsdb.pasteur.fr/…` link, sign in to your Pasteur
account and approve the application for the Klebsiella sequence-definition
database. Copy the verification code shown by Pasteur. Then complete the flow:

```bash
pixi run chronoclade authorize-pasteur \
  --destination ~/.local/share/chronoclade/typing --complete
```

Enter the verification code at the hidden prompt. There is no verifier command
argument, and the code is not printed or written to the access cache. Complete
the request within one hour; expired pending requests require a fresh start.
Your account must be registered for the requested database.

Build the references using the same destination and original key-file path:

```bash
pixi run chronoclade setup-typing --species klebsiella --prepare-db \
  --destination ~/.local/share/chronoclade/typing \
  --pasteur-secrets /path/to/pasteur-client-key.txt
```

Reference download/indexing is a separate step and can take substantial time.
Successful authorisation alone does not mean the references are installed or
that native typing is ready. Check `typing_setup.json` and `query_config.json`
after the build. See [query typing](query-typing.md) and
[input preparation](prepare.md).

## Private files and recovery

The destination's `credential-cache/` directory is owner-only (`0700`). Files
created by the authorisation flow are owner-readable/writable (`0600`):

- `pasteur-authorization-consumer.json`: the normalised consumer credentials.
- `pasteur-authorization-pending.json`: a temporary request, tied to that consumer,
  its creation time and a local nonce. Removed after a successful saved exchange.
- `plincer.json`: access tokens used by the upstream typing tools.

These are credentials, not analysis outputs. Keep them out of reports, shared
output bundles and version control. OAuth requests run in the installed plincer
environment; no additional core ChronoClade OAuth dependencies are needed.

Repeating the begin command reuses a valid pending request. If the browser says a
request has already been redeemed, the flow was interrupted after exchange, or
you need new consent, add `--restart` to the begin command. This obtains a fresh
request without deleting the existing access cache. Do not combine `--restart`
with `--complete`.

An unsuccessful exchange leaves existing access-cache contents intact. A
successful exchange replaces only the Klebsiella access grant and removes its
old session/request entries so the upstream tool obtains a fresh session. Other
hosts and databases remain preserved. Network responses and subprocess logs are
not displayed because they can contain credential material. Check account
registration, connectivity and the chosen installation destination before retrying.
