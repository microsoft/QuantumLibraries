// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using Microsoft.Jupyter.Core;
using Microsoft.Extensions.Logging;
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

        public EncodedData? Encode(object displayable) => displayable switch
        {
            CallSites sites => $@"Calls to {sites.Subject}:\n\n{
                    String.Join("\n---\n",
                        sites.Sites.Select(
                            call => String.Join("\n",
                                call.Select(
                                    frame => $"- {frame}"
                                )
                            )
                        )
                    )
                }".ToEncodedData(),
            _ => null
        };
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

        public EncodedData? Encode(object displayable) => displayable switch
        {
            CallSites sites => $@"
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
                ".ToEncodedData(),
            _ => null
        };

    }

}
