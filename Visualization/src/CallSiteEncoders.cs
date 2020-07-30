// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System.Text;
using Microsoft.Jupyter.Core;
using Microsoft.Extensions.Logging;
using Newtonsoft.Json;
using static NumSharp.Slice;
using Microsoft.Quantum.IQSharp.Jupyter;
using System.Linq;
using System;

namespace Microsoft.Quantum.Diagnostics.Emulation
{
    public class CallSitesToTextEncoder : IResultEncoder
    {
        private ILogger<CallSitesToTextEncoder> logger;
        public string MimeType => MimeTypes.PlainText;

        public CallSitesToTextEncoder(ILogger<CallSitesToTextEncoder> logger)
        {
            this.logger = logger;
        }

        public EncodedData? Encode(object displayable)
        {
            if (displayable is CallSites sites)
            {
                var list = String.Join("\n---\n",
                    sites.Sites.Select(
                        call => String.Join("\n",
                            call.Select(
                                frame => $"- {frame}"
                            )
                        )
                    )
                );
                return $"Calls to {sites.Subject}:\n\n{list}".ToEncodedData();
            }
            else return null;
        }

    }

    public class CallSitesToHtmlEncoder : IResultEncoder
    {
        private readonly ILogger<CallSitesToHtmlEncoder> logger;
        private readonly IConfigurationSource configurationSource;
        public string MimeType => MimeTypes.Html;

        public CallSitesToHtmlEncoder(
            ILogger<CallSitesToHtmlEncoder> logger,
            IConfigurationSource configurationSource
        )
        {
            this.logger = logger;
            this.configurationSource = configurationSource;
        }

        public EncodedData? Encode(object displayable)
        {
            if (displayable is CallSites sites)
            {
                logger.LogDebug("Using HTML encoder to output call site record.");
                return $@"
                    <details>
                        <summary>Calls to {sites.Subject}:</summary>
                        <ul>{
                            String.Join("\n",
                                sites.Sites.Select(
                                    call => $@"
                                        <li><ul>{
                                            String.Join("\n", call.Select(
                                                frame => $"<li>{frame}</li>"
                                            ))
                                        }</ul></li>
                                    "
                                )
                            )
                        }</ul>
                    </details>
                ".ToEncodedData();
            }
            else return null;
        }

    }

}
